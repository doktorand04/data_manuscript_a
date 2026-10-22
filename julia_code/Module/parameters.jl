# Parameter container and preallocated work arrays for the ADRD model.
#
# This file is meant to be `include`d from ADRD.jl, which should load
# PreallocationTools. `calc_spectral_D22` is defined elsewhere in the module.

"""
    EcoEvoParameters

Immutable container for all model parameters. Because it is immutable, "changing"
a parameter creates a new object, e.g. `p = @set p.S = p.S + 1`.

# Fields
- `s`, `a₀`: TODO describe
- `b`: TODO describe (default `4.0`)
- `α`, `β`, `γ`: TODO describe (defaults `0.5`, `0.1`, `1.0`)
- `δ`: diffusion coefficient of the consumers
- `δᵣ`: diffusion coefficient of the resources
- `N`: number of grid points per side (the spatial grid has `N^2` cells)
- `K₁ₛ`, `K₂ₛ`: inflow rates of resources 1 and 2, length `N^2`
- `δΔ`, `δᵣΔ`: consumer / resource diffusion operators (`δ` and `δᵣ` times the
  spectral Laplacian). Derived from `δ`, `δᵣ` and `N`; do not edit separately.
- `ess`: `true` once an evolutionarily stable state has been reached
- `ϵ_evo`: TODO describe (default `1e-6`)
- `S`: current number of consumer species

Use the positional constructor
`EcoEvoParameters(s, a₀, K₁ₛ, K₂ₛ, δ, δᵣ, α, β, γ)` for normal runs. Use the keyword
constructor to change `b`, `ϵ_evo` or other defaults.
"""
Base.@kwdef struct EcoEvoParameters
    s::Float64
    a₀::Float64
    b::Float64 = 4.0
    α::Float64 = 0.5
    β::Float64 = 0.1
    γ::Float64 = 1.0
    δ::Float64
    δᵣ::Float64
    N::Int
    K₁ₛ::Vector{Float64}
    K₂ₛ::Vector{Float64}
    δΔ::Matrix{Float64}  = δ  * calc_spectral_D22(N; L = 1.0)
    δᵣΔ::Matrix{Float64} = δᵣ * calc_spectral_D22(N; L = 1.0)
    ess::Bool = false
    ϵ_evo::Float64 = 1e-6
    S::Int = 1
end

# Positional constructor. The grid size `N` is inferred from the length of `K₁ₛ`, and
# the spectral Laplacian is computed once and shared by both diffusion operators.
function EcoEvoParameters(s, a₀, K₁ₛ, K₂ₛ, δ, δᵣ, α, β, γ)
    N = isqrt(length(K₁ₛ))
    N^2 == length(K₁ₛ) ||
        throw(ArgumentError("length(K₁ₛ) = $(length(K₁ₛ)) is not a perfect square"))
    length(K₂ₛ) == length(K₁ₛ) ||
        throw(DimensionMismatch("K₁ₛ and K₂ₛ must have the same length"))

    D22 = calc_spectral_D22(N; L = 1.0)
    return EcoEvoParameters(; s, a₀, α, β, γ, N, K₁ₛ, K₂ₛ, δ, δᵣ,
                            δΔ = δ * D22, δᵣΔ = δᵣ * D22)
end

"""
    ODECache

Preallocated work arrays for the diffusion terms, so the ODE right-hand side does
not allocate on every call. Each field is a `PreallocationTools.DiffCache`, which
also works with automatic differentiation (used by implicit solvers).
"""
struct ODECache{Tu,TR₁,TR₂}
    δΔuᵢ::Tu
    δᵣΔR₁::TR₁
    δᵣΔR₂::TR₂
end

function make_ode_cache(uₛ, R₁ₛ, R₂ₛ)
    return ODECache(
        PreallocationTools.DiffCache(uₛ),
        PreallocationTools.DiffCache(R₁ₛ),
        PreallocationTools.DiffCache(R₂ₛ),
    )
end