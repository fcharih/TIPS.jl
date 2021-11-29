function mutate!(p::Peptide, mutation_rate::Float)
    residues = split(p.sequence, "")
    for i in eachindex(residues)
        if rand(Float32) < mutation_rate
            residues[i] = rand(AMINO_ACIDS)
        end
    end
    p.sequence = join(residues)
end

function crossover!(p1::Peptide, p2::Peptide)
    s1 = split(p1.sequence, "")
    s2 = split(p2.sequence, "")
    new_s1, new_s2 = TPX(s1, s2)
    p1.sequence = join(new_s1)
    p2.sequence = join(new_s2)
end


