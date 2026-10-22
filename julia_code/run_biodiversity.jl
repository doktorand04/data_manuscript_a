# Generate biodiversity simulations using the ADRD module.
#
# Example usage in Julia REPL:
# empty!(ARGS)
# append!(ARGS, ["--a0","6.5","--s","1.0","--delta","0.0001","--nl","1","--nr","32","--alpha","0.5","--beta","0.1","--gamma","1.0","--output-dir","results"])
# include("run_biodiversity.jl")


using ArgParse
using OrdinaryDiffEq   # KenCarp4
using Sundials         # CVODE_BDF
using LinearSolve      # KrylovJL_GMRES

include(joinpath(@__DIR__, "Module", "ADRD.jl"))
using .ADRD

function parse_commandline()
    settings = ArgParseSettings(
        description = "Run the ADRD biodiversity generator.",
    )

    @add_arg_table! settings begin
        "--a0"
            help = "TODO: describe a₀"
            arg_type = Float64
            required = true
        "--s"
            help = "TODO: describe s"
            arg_type = Float64
            required = true
        "--delta"
            help = "TODO: describe δ"
            arg_type = Float64
            required = true
        "--nl"
            help = "TODO: describe Nₗ"
            arg_type = Int
            required = true
        "--nr"
            help = "TODO: describe Nᵣ"
            arg_type = Int
            required = true
        "--alpha"
            help = "TODO: describe α"
            arg_type = Float64
            required = true
        "--beta"
            help = "TODO: describe β"
            arg_type = Float64
            required = true
        "--gamma"
            help = "TODO: describe γ"
            arg_type = Float64
            required = true
        "--output-dir"
            help = "Directory where results are written"
            arg_type = String
            default = normpath(joinpath(@__DIR__, "..", "affinity_tradeoff_substitutability"))
    end

    return parse_args(settings)
end

function main()
    args = parse_commandline()

    a₀ = args["a0"]
    s  = args["s"]
    δ  = args["delta"]
    Nₗ = args["nl"]
    Nᵣ = args["nr"]
    α  = args["alpha"]
    β  = args["beta"]
    γ  = args["gamma"]
    path = abspath(expanduser(args["output-dir"]))

    mkpath(path)
    println("a₀ = $a₀, s = $s")

    ADRD.biodiversity_generator(
        Nᵣ, Nₗ, a₀, s, δ, α, β, γ,
        CVODE_BDF(linear_solver = :GMRES),
        KenCarp4(linsolve = KrylovJL_GMRES()),
        path,
    )
end

main()