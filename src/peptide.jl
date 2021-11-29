
AMINO_ACIDS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

mutable struct Peptide
    identifier::String
    sequence::String
    fitness::Float32
    Peptide(sequence::String) = new("hello", sequence, -Inf32)
end

Base.show(io::IO, peptide::Peptide) = print(io, "Peptide(\"$(peptide.sequence)\")")


"""
Creates a random peptide of length `size`
"""
random_peptide(size::Int64) = Peptide(join(rand(AMINO_ACIDS, size), ""))
