function save_peptides(filepath::String, peptides::Vector{Peptide})
    open(filepath, "w") do io
        write(io, "id,sequence,target_score,target_rank,max_off_target_score,max_off_target_rank,fitness\n")
        for (index, peptide) in enumerate(peptides)
            write(io, "$(peptide.id),$(peptide.sequence),$(peptide.target_score),$(peptide.target_rank),$(peptide.max_off_target_score),$(peptide.max_off_target_rank),$(peptide.fitness)\n")
        end
    end
end

function update_summary(filepath::String, generation::Int64, duration::AbstractFloat, peptides::Vector{Peptide})
    open(filepath, "a") do io
        if generation == 1
            write(io, "generation,duration,fittest_peptide_id,sequence,target_score,target_rank,max_off_target_score,max_off_target_rank,fitness\n")
        end
        fittest = peptides[argmax(map(p -> p.fitness, peptides))]
        write(io, "$(generation),$(duration),$(fittest.id),$(fittest.sequence),$(fittest.target_score),$(fittest.target_rank),$(fittest.max_off_target_score),$(fittest.max_off_target_rank),$(fittest.fitness)\n")
    end
end
