using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()
@everywhere using ArgParse
@everywhere using Random

@everywhere using TIPS
@everywhere using SPRINT
@everywhere using SPRINT.DataTypes
@everywhere using Evolutionary

function run_trial(run_index::Int64, proteins::Vector{Protein}, hsps::Set{HSP}, training_pairs::Vector{Vector{String}}, args)

    run_dir = "$(args["output"])/run$(run_index)"
    mkdir(run_dir)
    mkdir("$(run_dir)/peptides")

    if args["save_matrices"]
        matrices_dir = "$(run_dir)/score_matrices"
        mkdir(matrices_dir)
    end

    @info "Initializing a population of $(args["population_size"]) peptides..."
    Random.seed!(rand(1:100))
    peptides = [random_peptide("0-$(i)", args["peptide_length"]) for i in 1:args["population_size"]]

    @info "Scoring the ancestral peptide population..."
    score_matrix = score_peptides(proteins,
                          [create_protein(i, p.id, p.sequence) for (i, p) in enumerate(peptides)],
                            hsps,
                            training_pairs;
                            smer_sim_threshold = args["smer_similarity_threshold"],
                            hsp_sim_threshold = args["hsp_similarity_threshold"],
                            kmer_size = args["kmer_size"])

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
        score_matrix = score_peptides(proteins,
                          [create_protein(i, p.id, p.sequence) for (i, p) in enumerate(offspring)],
                            hsps,
                            training_pairs;
                            smer_sim_threshold = args["smer_similarity_threshold"],
                            hsp_sim_threshold = args["hsp_similarity_threshold"],
                            kmer_size = args["kmer_size"])

        @info "Computing the fitness of the new peptides..."
        compute_fitness!(offspring, args["target_id"], score_matrix, fitness_function)

        duration = time() - start_time
        @info "ITERATION $(generation) COMPLETED IN $(duration) seconds!"

        @info "Logging the results..."
        update_summary("$(run_dir)/summary.csv", generation, duration, offspring)
        save_peptides_csv("$(run_dir)/peptides/gen$(generation).csv", offspring)
        save_peptides_fasta("$(run_dir)/peptides/gen$(generation).fasta", offspring)

        if args["save_matrices"]
            FileIO.save_matrix("$(matrices_dir)/gen$(generation).mat")
        end

        peptides = [copy(peptide) for peptide in offspring]
        generation += 1

        # TODO check for convergence
        fittest = maximum(map(p -> p.fitness, peptides))
        if fittest > fittest_to_date
            fittest_to_date = fittest 
        end

        if fittest - fittest_to_date > args["convergence_variation"]
            iterations_since_change = 0
        else
            iterations_since_change += 1
        end

        if iterations_since_change ≥ args["convergence_iterations"]
            @info "======== CONVERGED ========"
            return
        end
    end
end

function main(args)

    @info "Loading the proteome sequences..."
    proteins = FileIO.load_sequences(args["fasta"])

    @info "Loading the training pairs..."
    training_pairs = FileIO.load_training_pairs(args["training_pairs"])

    @info "Loading the HSPs extracted from the proteome..."
    hsps = FileIO.load_hsps(args["hsps"], proteins)
    @info "Loaded $(length(hsps)) HSPs..."

    @info "Creating a directory structure for the evolution results..."
    if isdir(args["output"])
        @error "$(args["output"]) already exists and cannot be overwritten."
        return
    end

    mkdir(args["output"])

    for run_index in 1:args["n_runs"]
        @info "Starting run $(run_index)..."
        run_trial(run_index, proteins, hsps, training_pairs, args)
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
        "--convergence_variation"
            arg_type = Float64
            default = 3.0
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

    main(args)
end
