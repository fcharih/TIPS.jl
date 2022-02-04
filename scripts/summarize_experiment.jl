using ArgParse
using UnicodePlots
using CSV
using DataFrames
using NaturalSort

function summarize_run(run_dir::String)
    peptides_dir = joinpath(run_dir, "peptides")
    # Evolution of mean fitness
    peptide_files = filter(f -> occursin(".csv", f), readdir(peptides_dir))
    df = vcat([DataFrame(CSV.File(joinpath(peptides_dir, f))) for f in peptide_files]...)
    println(df)
end

function main(args)
    runs = sort(readdir(args["input"]), lt=natural)
    for run in runs
        run_path = joinpath(args["input"], run)
        summarize_run(run_path)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input", "-i"
            arg_type = String
            required = true
            help = "Path to the directory corresponding to an experiment."
    end

    args = parse_args(ARGS, s)

    main(args)
end

