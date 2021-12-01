module TIPS

using Evolutionary
using SPRINT.DataTypes: ScoreMatrix

export Peptide, random_peptide, AMINO_ACIDS, mutate!, crossover!, save_peptides_csv, save_peptides_fasta, compute_fitness!, update_summary

include("./peptide.jl")
include("./operators.jl")
include("./fileio.jl")
include("./fitness.jl")

end # module
