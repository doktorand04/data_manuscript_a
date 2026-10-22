# Eco-evolutionary simulation driver and result I/O for the ADRD module.
#
# This file is meant to be `include`d from ADRD.jl, which should load:
#   DataFrames, XLSX, a package providing `@set` (Accessors or Setfield),
#   and the DifferentialEquations pieces used below
#   (`CallbackSet`, `TerminateSteadyState`).

## -----------------------------------------------------------------------------
## Constants
## -----------------------------------------------------------------------------

const INITIAL_CONSUMER_DENSITY = 0.05  # initial density of every resident consumer
const TRAIT_RANGE_MARGIN       = 0.01  # invader traits are sampled on [-1+m, 1-m]
const TRAIT_RANGE_LENGTH       = 201   # number of points in the fitness landscape grid

# Eco-evo loop settings
const ECO_STEADY_STATE_TOL    = 1e-10  # tolerance for TerminateSteadyState (ecology)
const ECO_ODE_TOL             = 1e-6   # reltol/abstol for the ecological ODE solve
const ECO_END_TIME_PER_SPECIES = 10000.0
const EVO_INITIAL_TIME        = 1e10   # doubled after every evolutionary solve
const SELECTION_GRADIENT_TOL  = 1e-10  # evolutionary steady state criterion
const SAVE_EVERY_N_EVO_STEPS  = 5      # checkpoint frequency

## -----------------------------------------------------------------------------
## Entry points
## -----------------------------------------------------------------------------

"""
    biodiversity_generator(Nᵣ, Nₗ, a₀, s, δ, α, β, γ, ecosolver, evosolver, savepath;
                             invader_density = 1e-2)

Start an eco-evolutionary simulation from a single resident species.
The resource landscape is read from
`savepath/resource_landscape_<Nₗ>/ResourceLandscape_grid_size=<Nᵣ>.xlsx`
(sheet "Sheet1", columns `K1` and `K2`).
"""
function biodiversity_generator(Nᵣ, Nₗ, a₀, s, δ, α, β, γ,
                                  ecosolver, evosolver, savepath;
                                  invader_density = 1e-2)
    K₁, K₂ = load_resource_landscape(savepath, Nₗ, Nᵣ)
    println("Resources Loaded")

    p = EcoEvoParameters(s, a₀, K₁, K₂, δ, δ, α, β, γ)

    # Initial conditions
    uₛ  = INITIAL_CONSUMER_DENSITY * ones(Nᵣ^2)
    R₁ₛ = K₁ / 2
    R₂ₛ = K₂ / 2
    χₛ  = [0.0]

    cache = make_ode_cache(copy(uₛ), copy(R₁ₛ), copy(R₂ₛ))
    χᵣₙ   = trait_grid()
    uₘ    = invader_density * ones(Nᵣ^2)  # invaders start as small, flat densities

    eco_evo_simulation(Nₗ, χₛ, χᵣₙ, uₛ, uₘ, R₁ₛ, R₂ₛ, p, cache,
                       ecosolver, evosolver, savepath)
end

"""
    biodiversity_generator_from_mid(Nₗ, y, p, ecosolver, evosolver, savepath;
                                    invader_density = 1e-3)

Resume a simulation from a saved state vector `y`.
"""
function biodiversity_generator_from_mid(Nₗ, y, p, ecosolver, evosolver, savepath;
                                         invader_density = 1e-2)
    R₁ₛ, R₂ₛ, uₛ, χₛ = get_state_variables_eco_evo(y, p)

    u_cache = INITIAL_CONSUMER_DENSITY * ones(p.N^2)
    cache   = make_ode_cache(copy(u_cache), copy(R₁ₛ), copy(R₂ₛ))
    χᵣₙ     = trait_grid()
    uₘ      = invader_density * ones(p.N^2)

    eco_evo_simulation(Nₗ, χₛ, χᵣₙ, uₛ, uₘ, R₁ₛ, R₂ₛ, p, cache,
                       ecosolver, evosolver, savepath)
end

## -----------------------------------------------------------------------------
## Main loop
## -----------------------------------------------------------------------------

"""
Alternate between ecological relaxation and evolutionary dynamics, adding the
maximally fit invader after each evolutionary steady state, until no invader
can establish (an evolutionarily stable state).
"""
function eco_evo_simulation(Nₗ::Int, χₛ, χᵣₙ, uₛ, uₘ, R₁ₛ, R₂ₛ, p, cache,
                            ecosolver, evosolver, savepath)
    ecoodeopts = (;
        callback = CallbackSet(TerminateSteadyState(ECO_STEADY_STATE_TOL,
                                                    ECO_STEADY_STATE_TOL)),
        reltol = ECO_ODE_TOL,
        abstol = ECO_ODE_TOL,
        save_everystep = false,
    )

    lnum = 0
    while true
        lnum += 1

        ## -- Ecological dynamics until steady state
        yₑ = get_state_vector_eco(R₁ₛ, R₂ₛ, uₛ)
        χₛ, yₑ, p = extinction(χₛ, yₑ, p)

        t_end  = ECO_END_TIME_PER_SPECIES * length(χₛ)
        sol_eco = run_eco(yₑ, χₛ, (p, cache), t_end;
                          ode_solver = ecosolver, ode_opts = ecoodeopts)
        yₑ = sol_eco[end]
        println("Finished Ecol Modelling $lnum")

        ## -- Evolutionary dynamics until the selection gradient vanishes
        evoodeopts = (;
            reltol = 1e-6,
            abstol = [fill(1e-10, length(χₛ)); fill(1e-7, length(yₑ))],
            save_everystep = false,
        )

        yₑₑ        = [yₑ; χₛ]
        ∂λ∂χ       = calc_max_sel_grad(yₑₑ, p)
        total_time = EVO_INITIAL_TIME
        num        = 0
        while ∂λ∂χ >= SELECTION_GRADIENT_TOL
            sol_eco_evo = run_eco_evo(yₑₑ, (p, cache), total_time;
                                      ode_solver = evosolver, ode_opts = evoodeopts)
            yₑₑ        = sol_eco_evo[end]
            ∂λ∂χ       = calc_max_sel_grad(yₑₑ, p)
            total_time = 2 * total_time

            println("sg = ", ∂λ∂χ)
            flush(stdout)

            num += 1
            if mod(num, SAVE_EVERY_N_EVO_STEPS) == 0
                save_the_results(Nₗ, yₑₑ, p, savepath)
            end
        end
        println("Finished Eco-Evo Modelling $lnum")

        ## -- Look for an invader
        χₘ = calc_maximally_fit_invader_trait(χᵣₙ, yₑₑ, p)
        if isempty(χₘ)
            p = @set p.ess = true
            save_the_results_2(Nₗ, yₑₑ, p, savepath)
            println("Evolutionarily stable state reached")
            break
        end

        R₁ₛ, R₂ₛ, uₛᵦ, χₛᵦ = get_state_variables_eco_evo(yₑₑ, p)
        uₛ = [uₛᵦ uₘ]
        χₛ = [χₛᵦ; χₘ]
        p  = @set p.S = p.S + 1
    end
    return nothing
end

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------

trait_grid() = range(-1 + TRAIT_RANGE_MARGIN, 1 - TRAIT_RANGE_MARGIN;
                     length = TRAIT_RANGE_LENGTH)

resource_folder(Nₗ) = "resource_landscape_$Nₗ"

# Note: the exact format (including the missing space after "alpha = ...,") is kept
# so that files written by earlier versions of the code are still found.
function result_filename(p; suffix = "")
    (; a₀, s) = p
    return "a0 = $a₀, s = $s, alpha = $(p.α),beta = $(p.β), gamma = $(p.γ), " *
           "delta = $(p.δ)$suffix.xlsx"
end

function load_resource_landscape(savepath, Nₗ, Nᵣ)
    excelpath = joinpath(savepath, resource_folder(Nₗ),
                         "ResourceLandscape_grid_size=$Nᵣ.xlsx")
    isfile(excelpath) || error("Resource landscape file not found: $excelpath")
    df = DataFrames.DataFrame(XLSX.readtable(excelpath, "Sheet1"))
    return Float64.(df.K1), Float64.(df.K2)
end

"""
Write the state vector `y` and parameters `p` to an .xlsx file at `path`
(overwriting any existing file), one sheet per quantity.
"""
function write_state_xlsx(path, y, p)
    R₁, R₂, U, χ = get_state_variables_eco_evo(y, p)
    (; a₀, s) = p

    parameters = DataFrames.DataFrame(a0 = a₀, s = s, b = p.b, alpha = p.α,
                                      beta = p.β, gamma = p.γ, deltaC = p.δ,
                                      deltaR = p.δᵣ, ess = p.ess)
    inflow    = DataFrames.DataFrame(K1 = p.K₁ₛ, K2 = p.K₂ₛ)
    traits    = DataFrames.DataFrame(traits = χ)
    consumers = DataFrames.DataFrame(U, :auto)
    resources = DataFrames.DataFrame(R1 = R₁, R2 = R₂)

    mkpath(dirname(path))
    XLSX.writetable(path,
                    "Parameters"         => parameters,
                    "Inflow Rates"       => inflow,
                    "Traits"             => traits,
                    "Consumer Densities" => consumers,
                    "Resource Densities" => resources;
                    overwrite = true)
    return nothing
end

## -----------------------------------------------------------------------------
## Saving results
## -----------------------------------------------------------------------------

"Checkpoint a (possibly non-stable) state during the evolutionary dynamics."
function save_the_results(Nₗ::Int, y, p, savepath)
    path = joinpath(savepath, resource_folder(Nₗ), "Evolution", result_filename(p))
    return write_state_xlsx(path, y, p)
end

save_the_results(Nₗ, y, χ, p, savepath) = save_the_results(Nₗ, [y; χ], p, savepath)

"Save the evolutionarily stable state and delete the corresponding checkpoint."
function save_the_results_2(Nₗ::Int, y, p, savepath)
    evolution_dir = joinpath(savepath, resource_folder(Nₗ), "Evolution")
    filename = result_filename(p)
    write_state_xlsx(joinpath(evolution_dir, "Evolutionary Stable State", filename), y, p)

    checkpoint = joinpath(evolution_dir, filename)
    isfile(checkpoint) && rm(checkpoint)
    return nothing
end

"Save the state used to compute the fitness landscape of consumer `trait`."
function save_the_results_for_fitness_landscape(Nₗ::Int, y, p, savepath, trait)
    filename = result_filename(p; suffix = ", consumer = $trait")
    path = joinpath(savepath, resource_folder(Nₗ), "Evolution", filename)
    return write_state_xlsx(path, y, p)
end

"""
Append the trait and total-biomass time series of a PDE solution to an .xlsx file
under `savepath/resource_landscape_<Nₗ>/Trait Biomass/`.
"""
function save_trait_biomass(p, PDE_solution, init_time, savepath, Nₗ::Int)
    t = PDE_solution.t .+ init_time

    solution = Array(PDE_solution)
    χidxs = (2 + p.S) * p.N^2 + 1 : (2 + p.S) * p.N^2 + p.S
    χₛ    = solution[χidxs, :]'

    uidxs = 2 * p.N^2 + 1 : (2 + p.S) * p.N^2
    uₛ    = reshape(solution[uidxs, :], p.N^2, p.S, length(t))
    mean_density = dropdims(sum(uₛ, dims = 1), dims = 1)' ./ p.N^2

    (; a₀, s) = p
    parameters = DataFrames.DataFrame(a0 = a₀, s = s, b = p.b, alpha = p.α,
                                      beta = p.β, gamma = p.γ, deltaC = p.δ,
                                      deltaR = p.δᵣ, ess = p.ess)
    times     = DataFrames.DataFrame(time = Array(t))
    traits    = DataFrames.DataFrame(χₛ, :auto)
    densities = DataFrames.DataFrame(mean_density, :auto)

    folder = joinpath(savepath, resource_folder(Nₗ), "Trait Biomass")
    mkpath(folder)
    path = joinpath(folder, result_filename(p; suffix = ", consumer = $(p.S)"))

    # If the file already exists, append the new time series to the old one.
    # (This re-reads the whole file each call, so it gets slow for long runs.)
    if isfile(path)
        times     = append!(DataFrames.DataFrame(XLSX.readtable(path, "Time")), times)
        traits    = append!(DataFrames.DataFrame(XLSX.readtable(path, "Traits")), traits)
        densities = append!(DataFrames.DataFrame(XLSX.readtable(path, "Consumer Densities")), densities)
    end

    XLSX.writetable(path,
                    "Parameters"         => parameters,
                    "Time"               => times,
                    "Traits"             => traits,
                    "Consumer Densities" => densities;
                    overwrite = true)
    return nothing
end