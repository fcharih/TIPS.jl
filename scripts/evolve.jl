using ArgParse
using Random
using TIPS
using SPRINT

function main(args)

    @info "Loading the proteome sequences..."
    proteins = FileIO.load_sequences(args["protein_sequences"])

    @info "Loading the HSPs extracted from the proteome..."
    hsps = FileIO.load_hsps(args["hsps"])
    @info "Loaded $(length(hsps)) HSPs..."

    @info "Loading the training pairs..."
    training_pairs = FileIO.load_training_pairs(args["training_pairs"])

    @info "Initializing a population of $(args["population_size"]) peptides..."
    Random.seed!(args["seed"])
    peptides = [random_peptide(args["peptide_length"]) for i in 1:args["population_size"]]

    @info "Scoring the ancestral peptide population..."
    score_matrix = score_peptides(proteins,
                            [(p.identifier, p.sequence) for p in peptides],
                            hsps,
                            training_pairs;
                            smer_sim_threshold = args["smer_similarity_threshold"],
                            hsp_sim_threshold = args["hsp_similarity_threshold"],
                            kmer_size = args["kmer_size"])

    @info "Computing the fitness of ancestral peptides using computed interaction scores..."
    @Thread.threads for peptide in peptides
        peptide.fitness = compute_fitness(peptide.identifier, score_matrix, args["fitness"])
    end

    generation = 1
    while generation < args["max_generations"]
        @info "ITERATION $(generation) STARTED..."
        start_time = time()

        @info "Applying genetic operators to generate the next generation..."



        @info "Computing the interaction scores of newly generated peptides..."
        score_matrix = score_peptides(proteins,
                            [(p.identifier, p.sequence) for p in peptides],
                            hsps,
                            training_pairs;
                            smer_sim_threshold = args["smer_similarity_threshold"],
                            hsp_sim_threshold = args["hsp_similarity_threshold"],
                            kmer_size = args["kmer_size"])

        @info "Computing the fitness of the new peptides..."
        @Thread.threads for peptide in peptides
            peptide.fitness = compute_fitness(peptide.identifier, score_matrix, args["fitness_function"])
        end

        @info "Logging the results..."

        @info "ITERATION $(generation) COMPLETE!"
        generation += 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings()
    @add_arg_table s begin
        "--seed", "-s"
            arg_type = Int64
            required = false
            default = 1
            help = "Seed to initialize the population."
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
        "--output", "-o"
            arg_type = String
            required = true
            help = "Path to the directory where the experiment results should be saved."
    end

    args = parse_args(ARGS, s)

    main(args)
end
