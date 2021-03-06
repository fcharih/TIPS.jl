using FASTX

function save_peptides_csv(filepath::String, peptides::Vector{Peptide})
    # Write a CSV file
    open(filepath, "w") do io
        write(io, "id,sequence,target_score,target_rank,max_off_target_score,fitness\n")
        for (index, peptide) in enumerate(peptides)
            write(io, "$(peptide.id),$(peptide.sequence),$(peptide.target_score),$(peptide.target_rank),$(peptide.max_off_target_score),$(peptide.fitness)\n")
        end
    end
end

function save_peptides_fasta(filepath::String, peptides::Vector{Peptide})
    sorted = sort(peptides, by = p -> p.id)
    open(filepath, "w") do io
        for peptide in peptides
            write(io, ">$(peptide.id)\n$(peptide.sequence)\n")
        end
    end
end

function update_summary(filepath::String, generation::Int64, duration::AbstractFloat, peptides::Vector{Peptide}, seed::Integer)
    open(filepath, "a") do io
        if generation == 1
            write(io, "generation,seed,duration,fittest_peptide_id,sequence,target_score,target_rank,max_off_target_score,fitness\n")
        end
        fittest = peptides[argmax(map(p -> p.fitness, peptides))]
        write(io, "$(generation),$(seed),$(duration),$(fittest.id),$(fittest.sequence),$(fittest.target_score),$(fittest.target_rank),$(fittest.max_off_target_score),$(fittest.fitness)\n")
    end
end
