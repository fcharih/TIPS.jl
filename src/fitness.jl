using StatsBase: competerank

function compute_fitness!(peptides::Vector{Peptide}, target_id::String, score_matrix::ScoreMatrix, fitness_function::Function)

    @Threads.threads for peptide in peptides
        target_index = score_matrix.protein_ids[target_id]
        peptide_index = score_matrix.protein_ids[peptide.id]
        
        # Compute the metrics
        target_score = score_matrix.matrix[peptide_index, target_index]
        target_rank = competerank(score_matrix.matrix[peptide_index, :], rev = true)[target_index]
        max_off_target_score, max_off_target_index = findmax(score_matrix.matrix[peptide_index, 1:end .!= target_index])

        # Need to account for the fact that the index in the slice that excludes the target may shift the
        # position of max_off_target by -1 if max_off_target_index ≥ target_index
        max_off_target_index_with_target = max_off_target_index < target_index ? max_off_target_index : max_off_target_index + 1
        max_off_target_rank = competerank(score_matrix.matrix[peptide_index, :], rev = true)[max_off_target_index_with_target]

        # Set the peptide data
        peptide.target_score = target_score
        peptide.target_rank = target_rank
        peptide.max_off_target_score = max_off_target_score
        peptide.max_off_target_rank = max_off_target_rank

        # Compute the fitness
        # (Abreviations are useful when defining the fitness function in the command line call)
        n_prot = size(score_matrix.matrix)[1] - length(peptides)
        ts = target_score
        tr = target_rank
        mots = max_off_target_score
        motr = max_off_target_rank
        peptide.fitness = fitness_function(ts, tr, mots, motr, n_prot)
    end
end
