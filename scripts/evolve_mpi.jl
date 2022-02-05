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

struct SequencesAndHsps
    sequences::Vector{Protein}
    hsps::Set{HSP}
end

function dispatch_hsp_extraction(peptides, comm)
    for slave_id in 1:MPI.Comm_size(comm)-1
        MPI.send(peptides, slave_id, MPI.MPI_ANY_TAG, comm)
    end

    hsps = Extraction.extract_hsps(
        proteins;
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

    for slave_id in 1:MPI.Comm_size(comm)-1
        peptide_hsps, status = MPI.recv(slave_id, MPI.MPI_ANY_TAG, comm)
        union!(hsps, peptide_hsps)
    end

    hsps
end

function dispatch_prediction(proteins, peptides, hsps, training_pairs, comm, args)
    
    payload = SequencesAndHsps([proteins..., peptides...], hsps)

    for slave_id in 1:MPI.Comm_size(comm)-1
        MPI.send(paylod, slave_id, PREDICTION_TAG, comm)
    end

    matrix = Prediction.score_interactions(
        payload.sequences,
        hsps,
        training_pairs[1:MPI.Comm_size(comm):end],
        args["kmer_size"]
    )

    for slave_id in 1:MPI.Comm_size(comm)-1
        sub_matrix, status = MPI.recv(slave_id, MPI.MPI_ANY_TAG, comm)
        matrix += sub_matrix
    end

    matrix
end


function master_process(comm, args)

    rank = MPI.Comm_rank(comm)

    @info "Loading the proteome sequences..."
    proteins = FileIO.load_sequences(args["fasta"])

    @info "Loading the training pairs..."
    training_pairs = FileIO.load_training_pairs(args["training_pairs"])

    @info "Loading the HSPs extracted from the proteome..."
    training_hsps = FileIO.load_hsps(args["hsps"], proteins)
    @info "Loaded $(length(hsps)) HSPs..."

    @info "Creating a directory structure for the evolution results..."
    if isdir(args["output"])
        @error "$(args["output"]) already exists and cannot be overwritten."
        return
    end

    mkdir(args["output"])

    run_dir = "$(args["output"])/run$(run_index)"
    summary_file = "$(run_dir)/summary.csv"
    mkdir(run_dir)
    mkdir("$(run_dir)/peptides")

    if args["save_matrices"]
        matrices_dir = "$(run_dir)/score_matrices"
        mkdir(matrices_dir)
    end

    @info "Initializing a population of $(args["population_size"]) peptides..."
    seed = rand(1:100)
    Random.seed!(seed)
    peptides = [random_peptide("0-$(i)", args["peptide_length"]) for i in 1:args["population_size"]]

    hsps = dispatch_hsp_extraction(peptides, comm)
    score_matrix = dispatch_prediction(proteins, peptides, hsps, training_pairs, comm, args)

end

function slave_process(comm, args)
    rank = MPI.Comm_rank(comm)
    world_size = MPI.Comm_size(comm)

    @info "Loading the proteome sequences..."
    proteins = FileIO.load_sequences(args["fasta"])

    @info "Loading the training pairs..."
    training_pairs = FileIO.load_training_pairs(args["training_pairs"])

    @info "Loading the HSPs extracted from the proteome..."
    training_hsps = FileIO.load_hsps(args["hsps"], proteins)
    @info "Loaded $(length(hsps)) HSPs..."

    while true
        payload, status = MPI.recv(0, MPI.MPI_ANY_TAG, comm)

        if status.tag == EXTRACTION_TAG
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
            MPI.send(hsps, 0, 1, comm)
        else
            matrix = Prediction.score_interactions(
                payload.sequences,
                payload.hsps,
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
        master_process(comm, args)
    else
        slave_process(comm, args)
    end
end
