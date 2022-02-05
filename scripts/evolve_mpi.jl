using ArgParse
using MPI
using Random
using TIPS

using Evolutionary
using CSV
using DataFrames

using SPRINT
using SPRINT.DataTypes

EXTRACTION_TAG = 1001
PREDICTION_TAG = 1002

struct Payload
    sequences::Vector{Protein}
    hsps::Set{HSP}
end

function dispatch_hsp_extraction(sequences, comm)
    for worker_id in 1:MPI.Comm_size(comm)-1
        MPI.send(sequences, worker_id, EXTRACTION_TAG, comm)
        @info "Dispatched sequences to worker $(worker_id)..."
    end

    hsps = Extraction.extract_hsps(
        sequences;
        smer_sim_threshold = args["smer_similarity_threshold"],
        hsp_sim_threshold = args["hsp_similarity_threshold"],
        kmer_size = args["kmer_size"],
        new_sequences_only = true,
        trivial_hsps = false,
        verbose = false,
        process_id = 1,
        num_processes = MPI.Comm_size(comm),
        as_matrix = false
    )

    @info "Main process done extracting hsps."

    for worker_id in 1:MPI.Comm_size(comm)-1
        peptide_hsps, status = MPI.recv(worker_id, MPI.MPI_ANY_TAG, comm)
        @info "Main processed received HSPs from worker $(worker_id)..."
        union!(hsps, peptide_hsps)
    end

    hsps
end

function dispatch_prediction(sequences, training_hsps, peptide_hsps, training_pairs, comm, args)
    
    payload = Payload(sequences, peptide_hsps)

    for worker_id in 1:MPI.Comm_size(comm)-1
        MPI.send(payload, worker_id, PREDICTION_TAG, comm)
    end

    matrix = Prediction.score_interactions(
        payload.sequences,
        union(peptide_hsps, training_hsps),
        training_pairs[1:MPI.Comm_size(comm):end],
        args["kmer_size"]
    )

    for worker_id in 1:MPI.Comm_size(comm)-1
        sub_matrix, status = MPI.recv(worker_id, MPI.MPI_ANY_TAG, comm)
        matrix += sub_matrix
    end

    protein_ids = Dict(p.protein_id => i for (i, p) in enumerate(sequences))

    ScoreMatrix(protein_ids, matrix)
end

function isconverged(summary_file, args)::Bool
    evolution = CSV.File(summary_file) |> DataFrame

    if nrow(evolution) < args["min_generations"]
        return false
    end

    # Compute fitness delta
    Δfitness = []
    for i in 2:nrow(evolution)
        push!(Δfitness, evolution.fitness[i] - evolution.fitness[i-1])
    end

    if all(x -> x < args["convergence_variation"], Δfitness[end - args["convergence_iterations"]:end])
        return true
    end
end

function run(run_id, comm, args)

    run_dir = "$(output_dir)/run$(run_id)"
    summary_file = "$(run_dir)/summary.csv"
    mkdir("$(run_dir)/peptides")

    if args["save_matrices"]
        matrices_dir = "$(run_dir)/score_matrices"
        mkdir(matrices_dir)
    end

    @info "Initializing a population of $(args["population_size"]) peptides..."
    seed = rand(1:100)
    Random.seed!(seed)
    peptides = [random_peptide("0-$(i)", args["peptide_length"]) for i in 1:args["population_size"]]

    sequences = [deepcopy(proteins); [create_protein(p.id, p.sequence; new_sequence=true) for p in peptides]]
    peptide_hsps = dispatch_hsp_extraction(sequences, comm)
    score_matrix = dispatch_prediction(sequences, training_hsps, peptide_hsps, training_pairs, comm, args)

    fitness_expression = Meta.parse(args["fitness_function"])
    fitness_function = @eval function(ts, tr, mots, motr, n_prot)
        $fitness_expression
    end

    @info "Computing the fitness of ancestral peptides using computed interaction scores..."
    compute_fitness!(peptides, args["target_id"], score_matrix, fitness_function)

    # Initialize tournament
    tournament_function = tournament(args["tournament_size"], false)

    generation = 1
    fittest_to_date = -Inf32
    iterations_since_change = 0
    while true
        @info "ITERATION $(generation) STARTED..."
        start_time = time()

        @info "Applying genetic operators to generate the next generation..."
        offspring::Vector{Peptide} = []

        # Peptide selection
        fitnesses = [peptide.fitness for peptide in peptides]
        peptide_selection = tournament_function(fitnesses, args["population_size"])

        for (i, index_of_selected) in enumerate(peptide_selection)
            peptide_id = "$(generation)-$(length(offspring) + 1)"
            push!(offspring, Peptide(peptide_id, peptides[index_of_selected].sequence))
        end

        # Crossover
        shuffle!(offspring)
        for i in 1:2:floor(Int(length(offspring)/2))
            crossover!(offspring[i], offspring[i + 1])
        end

        # Mutate
        for peptide in offspring
            if rand(Float64) < args["mutation_probability"]
                mutate!(peptide, args["mutation_rate"])
            end
        end

        @info "Computing the interaction scores of newly generated peptides..."
        sequences = [deepcopy(proteins); [create_protein(p.id, p.sequence; new_sequence=true) for p in offspring]]
        peptide_hsps = dispatch_hsp_extraction(sequences, comm)
        score_matrix = dispatch_prediction(sequences, training_hsps, peptide_hsps, training_pairs, comm, args)

        @info "Computing the fitness of the new peptides..."
        compute_fitness!(offspring, args["target_id"], score_matrix, fitness_function)

        duration = time() - start_time
        @info "ITERATION $(generation) COMPLETED IN $(duration) seconds!"

        @info "Logging the results..."
        update_summary(summary_file, generation, duration, offspring, seed)
        save_peptides_csv("$(run_dir)/peptides/gen$(generation).csv", offspring)
        save_peptides_fasta("$(run_dir)/peptides/gen$(generation).fasta", offspring)

        if args["save_matrices"]
            FileIO.save_matrix("$(run_dir)/matrices/gen$(generation).mat")
        end

        peptides = [copy(peptide) for peptide in offspring]
        generation += 1

        # Check for convergence
        if isconverged(summary_file, args)
            @info "======== CONVERGED after $(generation - 1) generations ========"
            return
        end
    end
end

function main_process(comm, args)

    rank = MPI.Comm_rank(comm)

    @info "Loading the proteome sequences..."
    proteins = FileIO.load_sequences(args["fasta"])

    @info "Loading the training pairs..."
    training_pairs = FileIO.load_training_pairs(args["training_pairs"])

    @info "Loading the HSPs extracted from the proteome..."
    training_hsps = FileIO.load_hsps(args["hsps"], proteins)
    @info "Loaded $(length(training_hsps)) HSPs..."

    @info "Creating a directory structure for the evolution results..."
    if isdir(args["output"])
        @error "$(args["output"]) already exists and cannot be overwritten."
        return
    end

    output_dir = args["output"]
    mkdir(output_dir)

    for i in 1:args["n_runs"]
        run(i, comm, args)
    end
end

function worker_process(comm, args)
    rank = MPI.Comm_rank(comm)
    world_size = MPI.Comm_size(comm)

    @info "Loading the proteome sequences..."
    proteins = FileIO.load_sequences(args["fasta"])

    @info "Loading the training pairs..."
    training_pairs = FileIO.load_training_pairs(args["training_pairs"])

    @info "Loading the HSPs extracted from the proteome..."
    training_hsps = FileIO.load_hsps(args["hsps"], proteins)
    @info "Loaded $(length(training_hsps)) HSPs..."

    while true
        payload, status = MPI.recv(0, MPI.MPI_ANY_TAG, comm)

        if status.tag == EXTRACTION_TAG
            @info "Worker $(MPI.Comm_rank(comm)) received sequences for HSP extraction..."
            hsps = Extraction.extract_hsps(
                payload;
                smer_sim_threshold = args["smer_similarity_threshold"],
                hsp_sim_threshold = args["hsp_similarity_threshold"],
                kmer_size = args["kmer_size"],
                new_sequences_only = true,
                trivial_hsps = false,
                verbose = false,
                process_id = MPI.Comm_rank(comm) + 1,
                num_processes = MPI.Comm_size(comm),
                as_matrix = false
            )
            @info "Worker $(MPI.Comm_rank(comm)) sending HSPs back..."
            MPI.send(hsps, 0, 1, comm)
        else
            @info "Worker $(MPI.Comm_rank(comm)) received HSPs for scoring..."
            merged_hsps = union(payload.hsps, training_hsps)
            matrix = Prediction.score_interactions(
                payload.sequences,
                merged_hsps,
                training_pairs[rank+1:MPI.Comm_size(comm):end],
                args["kmer_size"]
            )
            MPI.send(matrix, 0, 1, comm)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings()
    @add_arg_table s begin
        "--n_runs", "-n"
            arg_type = Int64
            required = false
            default = 3
            help = "Number of initializations/runs..."
        "--target_id", "-i"
            arg_type = String
            required = true
            help = "ID of the target against which a peptide should be evolved."
        "--fasta", "-f"
            arg_type = String
            required = true
            help = "Path to the files containing the proteome sequences."
        "--hsps", "-s"
            arg_type = String
            required = true
            help = "Path to the files containing the processed HSPs."
        "--training_pairs", "-t"
            arg_type = String
            required = true
            help = "Path to the files containing the validated protein interactions."
        "--peptide_length", "-l"
            arg_type = Int64
            default = 20
            required = false
            help = "Length of peptide to evolve."
        "--population_size", "-p"
            arg_type = Int64
            default = 1000
            required = false
            help = "Population size (i.e. number of peptides to evolve)."
        "--min_generations"
            arg_type = Int64
            default = 50
            required = false
            help = "Minimum number of generations"
        "--convergence_variation"
            arg_type = Float64
            default = 1.0
            required = false
            help = "Width of score braket among which a run is considered to have converged after `convergence_iterations` iterations within that bracket."
        "--convergence_iterations"
            arg_type = Int64
            default = 10
            required = false
            help = "Number of iterations after which the run has converged if the score of fittest peptides remains bounded by `convergence_variation`."
        "--save_matrices"
            action = :store_true
            help = "Whether the score matrices for an iteration should be saved."
        "--kmer_size"
            arg_type = Int64
            default = 20
            required = false
            help = "Minimum HSP length."
        "--smer_similarity_threshold"
            arg_type = Int64
            default = 15
            required = false
            help = "Minimum threshold needed for 2 smers to be considered similar."
        "--hsp_similarity_threshold"
            arg_type = Int64
            default = 35
            required = false
            help = "Minimum threshold needed for 2 regions to be considered HSPs."
        "--crossover_probability"
            arg_type = Float64
            default = 0.6
            required = false
            help = "Probability of two peptides undergoing cross-over."
        "--mutation_probability"
            arg_type = Float64
            default = 0.5
            required = false
            help = "Probability of two peptides undergoing mutation."
        "--mutation_rate"
            arg_type = Float64
            default = 0.2
            required = false
            help = "Probability that a given residue in a mutating peptide will mutate."
        "--tournament_size"
            arg_type = Int64
            default = 3
            required = false
            help = "Number of peptides involved in a selection tournament."
        "--fitness_function"
            arg_type = String
            default = "mots > 140 ? 0 : ts * (1 - tr/n_prot)"
            required = false
            help = "Fitness function."
        "--output", "-o"
            arg_type = String
            required = true
            help = "Path to the directory where the experiment results should be saved."
    end

    args = parse_args(ARGS, s)

    MPI.Init()
    comm = MPI.COMM_WORLD

    if MPI.Comm_rank(comm) == 0
        main_process(comm, args)
    else
        worker_process(comm, args)
    end
end
