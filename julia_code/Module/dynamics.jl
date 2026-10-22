# ODE right-hand sides, solver wrappers and state-vector helpers.
#
# State vector layout (N = p.N grid points per side, S = p.S species):
#
#   y = [R₁ (N²); R₂ (N²); u₁ (N²); …; u_S (N²); χ₁, …, χ_S]
#
# The ecological system uses only the first (2 + S)·N² entries. The eco-evolutionary
# system appends the S traits χ.
#
# This file is meant to be `include`d from ADRD.jl, which should load
# OrdinaryDiffEq, SciMLBase, Sundials, LinearSolve, LinearAlgebra and
# PreallocationTools. `TRAIT_BOUNDARY_TOL` is defined in fitness.jl.

# Densities below this value are not allowed to decrease any further.
const ODE_EXTINCTION_LIMIT = 1e-16

## -----------------------------------------------------------------------------
## Solver wrappers
## -----------------------------------------------------------------------------

"""
    run_eco(y₀, χₛ, psys, tₑ; ode_solver, ode_opts)

Integrate the ecological dynamics from `y₀` over `(0, tₑ)` for fixed traits `χₛ`.
`psys = (p, cache)` holds the parameters and the preallocated work arrays.
Returns the ODE solution.
"""
function run_eco(y₀, χₛ, psys, tₑ;
                 ode_solver = CVODE_BDF(linear_solver = :GMRES),
                 ode_opts = (; reltol = 1e-6, abstol = 1e-6, save_everystep = false))
    ode_prob = ODEProblem(ode_sys_eco!, y₀, (0.0, tₑ), (psys, χₛ))
    return solve(ode_prob, ode_solver; ode_opts...)
end

"""
    run_eco_evo(y₀, psys, tₑ; ode_solver, ode_opts)

Integrate the coupled ecological and trait dynamics from `y₀` (which includes the
traits) over `(0, tₑ)`. `psys = (p, cache)`. Returns the ODE solution.
"""
function run_eco_evo(y₀, psys, tₑ;
                     ode_solver = KenCarp4(linsolve = KrylovJL_GMRES()),
                     ode_opts = (; reltol = 1e-6, abstol = 1e-6, save_everystep = true))
    ode_prob = ODEProblem{true, SciMLBase.AutoSpecialize}(ode_sys_eco_evo!, y₀,
                                                          (0.0, tₑ), psys)
    return solve(ode_prob, ode_solver; ode_opts...)
end

## -----------------------------------------------------------------------------
## Right-hand sides
## -----------------------------------------------------------------------------

"""
    ode_sys_eco!(dy, y, (psys, χₛ), t)

In-place right-hand side of the ecological dynamics for fixed traits `χₛ`.

Resources: `dRᵢ/dt = γ(Kᵢ - Rᵢ) + δᵣΔRᵢ - Σₛ Cᵢ G uₛ`.
Consumers: `duₛ/dt = (G - μ) uₛ + δΔuₛ`, with background mortality `μ = 1/α - 1`.
Derivatives are clipped at zero for densities below `ODE_EXTINCTION_LIMIT`, so
densities cannot become negative.
"""
function ode_sys_eco!(dy, y, (psys, χₛ), t)
    p, cache = psys

    R₁ₛ, R₂ₛ, uₛ = get_state_variables_eco(y, p)
    dR₁ₛ, dR₂ₛ, duₛ = get_state_variables_eco(dy, p)

    # Background mortality: chosen so that the net growth rate G - μ approaches 1
    # as T → ∞ (G → 1/α).
    μ = 1 / p.α - 1

    # Resource diffusion
    δᵣΔR₁ = PreallocationTools.get_tmp(cache.δᵣΔR₁, R₁ₛ)
    δᵣΔR₂ = PreallocationTools.get_tmp(cache.δᵣΔR₂, R₂ₛ)
    mul!(δᵣΔR₁, p.δᵣΔ, R₁ₛ)
    mul!(δᵣΔR₂, p.δᵣΔ, R₂ₛ)

    # Resource renewal and diffusion
    for j ∈ eachindex(R₁ₛ, R₂ₛ, dR₁ₛ, dR₂ₛ, p.K₁ₛ, p.K₂ₛ)
        dR₁ₛ[j] = p.γ * (p.K₁ₛ[j] - R₁ₛ[j]) + δᵣΔR₁[j]
        dR₂ₛ[j] = p.γ * (p.K₂ₛ[j] - R₂ₛ[j]) + δᵣΔR₂[j]
    end

    # Consumer dynamics and the effect of consumption on the resources
    for i ∈ axes(uₛ, 2)
        uᵢ  = view(uₛ, :, i)
        duᵢ = view(duₛ, :, i)

        # Consumer diffusion
        δΔuᵢ = PreallocationTools.get_tmp(cache.δΔuᵢ, uᵢ)
        mul!(δΔuᵢ, p.δΔ, uᵢ)

        # Affinities are constant across space, so compute them once per species
        a₁, a₂ = calc_a1_a2(χₛ[i], p)

        for j ∈ eachindex(uᵢ, duᵢ, R₁ₛ, dR₁ₛ, R₂ₛ, dR₂ₛ)
            X₁ = max(0.0, a₁ * R₁ₛ[j])
            X₂ = max(0.0, a₂ * R₂ₛ[j])

            C₁, C₂ = calc_C1_C2(X₁, X₂, p)
            G = calc_G(X₁, X₂, p)

            duᵢ[j]   = (G - μ) * uᵢ[j] + δΔuᵢ[j]
            dR₁ₛ[j] -= C₁ * G * uᵢ[j]
            dR₂ₛ[j] -= C₂ * G * uᵢ[j]

            duᵢ[j] = uᵢ[j] >= ODE_EXTINCTION_LIMIT ? duᵢ[j] : max(0.0, duᵢ[j])
        end
    end

    # Keep resource densities non-negative
    for j ∈ eachindex(R₁ₛ, R₂ₛ, dR₁ₛ, dR₂ₛ)
        dR₁ₛ[j] = R₁ₛ[j] >= ODE_EXTINCTION_LIMIT ? dR₁ₛ[j] : max(0.0, dR₁ₛ[j])
        dR₂ₛ[j] = R₂ₛ[j] >= ODE_EXTINCTION_LIMIT ? dR₂ₛ[j] : max(0.0, dR₂ₛ[j])
    end

    return nothing
end

"""
    ode_sys_eco_evo!(dy, y, psys, t)

In-place right-hand side of the coupled eco-evolutionary dynamics. The ecological
part is [`ode_sys_eco!`](@ref). Each trait follows `dχₛ/dt = ϵ_evo · ∂λ/∂χ`, where
the selection gradient is the `u²`-weighted average of the local gradients of `G`.
A trait at the boundary ±1 is not allowed to move outwards.
"""
function ode_sys_eco_evo!(dy, y, psys, t)
    p, cache = psys

    R₁ₛ, R₂ₛ, uₛ, χₛ = get_state_variables_eco_evo(y, p)
    _, _, _, dχₛ = get_state_variables_eco_evo(dy, p)

    # Ecological dynamics (first (2 + S)·N² entries)
    y_eco  = view(y, 1:(2 + p.S) * p.N^2)
    dy_eco = view(dy, 1:(2 + p.S) * p.N^2)
    ode_sys_eco!(dy_eco, y_eco, (psys, χₛ), t)

    # Selection gradient dynamics
    for i ∈ eachindex(χₛ, dχₛ)
        uᵢ = view(uₛ, :, i)
        ∂λ∂χ = sum((uᵢ[j]^2) * calc_local_sel_grad(χₛ[i], max(0.0, R₁ₛ[j]),
                                                    max(0.0, R₂ₛ[j]), p)
                   for j ∈ eachindex(R₁ₛ, R₂ₛ, uᵢ))
        ∂λ∂χ /= sum(uᵢ[j]^2 for j ∈ eachindex(uᵢ))

        if (1 - abs(χₛ[i])) >= TRAIT_BOUNDARY_TOL
            dχₛ[i] = p.ϵ_evo * ∂λ∂χ
        elseif χₛ[i] > 1.0 - TRAIT_BOUNDARY_TOL
            dχₛ[i] = min(0.0, p.ϵ_evo * ∂λ∂χ)
        elseif χₛ[i] < -(1.0 - TRAIT_BOUNDARY_TOL)
            dχₛ[i] = max(0.0, p.ϵ_evo * ∂λ∂χ)
        end
    end

    return nothing
end

## -----------------------------------------------------------------------------
## State vector helpers (see the layout at the top of this file)
## -----------------------------------------------------------------------------

"Return views `R₁ₛ, R₂ₛ, uₛ` into the ecological state vector `y` (`uₛ` is `N² × S`)."
function get_state_variables_eco(y, p)
    R₁ₛ = view(y, 1:p.N^2)
    R₂ₛ = view(y, p.N^2+1:2*p.N^2)

    uidxs = 2*p.N^2+1:(2+p.S)*p.N^2
    uₛ = reshape(view(y, uidxs), p.N^2, p.S)

    return R₁ₛ, R₂ₛ, uₛ
end

"Return views `R₁ₛ, R₂ₛ, uₛ, χₛ` into the eco-evolutionary state vector `y`."
function get_state_variables_eco_evo(y, p)
    R₁ₛ = view(y, 1:p.N^2)
    R₂ₛ = view(y, p.N^2+1:2*p.N^2)

    uidxs = 2*p.N^2+1:(2+p.S)*p.N^2
    uₛ = reshape(view(y, uidxs), p.N^2, p.S)

    χidxs = (2+p.S)*p.N^2+1:(2+p.S)*p.N^2+p.S
    χₛ = view(y, χidxs)

    return R₁ₛ, R₂ₛ, uₛ, χₛ
end

"Assemble the ecological state vector from `R₁ₛ`, `R₂ₛ` and `uₛ`."
get_state_vector_eco(R₁ₛ, R₂ₛ, uₛ) = [R₁ₛ; R₂ₛ; vec(uₛ)]

"Assemble the eco-evolutionary state vector from `R₁ₛ`, `R₂ₛ`, `uₛ` and `χₛ`."
get_state_vector_eco_evo(R₁ₛ, R₂ₛ, uₛ, χₛ) = [R₁ₛ; R₂ₛ; vec(uₛ); χₛ]