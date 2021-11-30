AMINO_ACIDS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

mutable struct Peptide
    id::String
    sequence::String
    target_score::Union{Float32, Nothing}
    target_rank::Union{Int64, Nothing}
    max_off_target_score::Union{Float32, Nothing}
    max_off_target_rank::Union{Int64, Nothing}
    fitness::Union{Float32, Nothing}
    Peptide(id::String, sequence::String, target_score::Union{Float32, Nothing}, target_rank::Union{Int64, Nothing}, max_off_target_score::Union{Float32, Nothing}, max_off_target_rank::Union{Int64, Nothing}, fitness::Union{Float32, Nothing}) = new(id, sequence, target_score, target_rank, max_off_target_score, max_off_target_rank, fitness)
    Peptide(id::String, sequence::String) = new(id, sequence, nothing, nothing, nothing, nothing, nothing)
end

Base.copy(p::Peptide) = Peptide(p.id, p.sequence, p.target_score, p.target_rank, p.max_off_target_score, p.max_off_target_rank, p.fitness)

Base.show(io::IO, peptide::Peptide) = print(io, "Peptide(\"$(peptide.sequence)\")")

"""
Creates a random peptide of length `size`
"""
random_peptide(id::String, size::Int64) = Peptide(id, join(rand(AMINO_ACIDS, size), ""))

"""
Checks if the peptide has been scored/i.e. fitness calculated.
"""
isscored(p::Peptide) = p.fitness != -Inf32
