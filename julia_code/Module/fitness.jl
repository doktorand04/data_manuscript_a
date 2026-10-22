# Trait-dependent resource affinities, the growth function G and selection gradients.
#
# This file is meant to be `include`d from ADRD.jl, which should load ForwardDiff.
# `get_state_variables_eco_evo` is defined elsewhere in the module.

# A trait χ closer than this to the boundary ±1 is treated as sitting on the boundary.
const TRAIT_BOUNDARY_TOL = 1e-6

## -----------------------------------------------------------------------------
## Affinities and growth function
## -----------------------------------------------------------------------------

"""
    calc_a1_a2(χ, p) -> (a₁, a₂)

Affinities for resources 1 and 2 of a consumer with trait `χ ∈ [-1, 1]`.

For `p.a₀ == 0` the trade-off is linear: `a₁ = 1 - χ`, `a₂ = 1 + χ`. Otherwise it is
nonlinear, with shape set by `p.a₀` and `p.b`. In both cases `a₁(0) = a₂(0) = 1`
and `a₁(1) = a₂(-1) = 0`.
"""
function calc_a1_a2(χ, p)
    (; a₀, b) = p
    if a₀ == 0
        a₁ = 1 - χ
        a₂ = 1 + χ
    else
        den = a₀ - b * log(0.5 * (exp(-a₀) + 1))
        a₁ₙ = b * log(0.5 * (exp(-a₀ * χ) + 1)) - a₀ * χ
        a₂ₙ = b * log(0.5 * (exp(a₀ * χ) + 1)) + a₀ * χ
        a₁ = 1 + a₁ₙ / den
        a₂ = 1 + a₂ₙ / den
    end
    return a₁, a₂
end

"""
    calc_T(X₁, X₂, p)
    calc_T(a₁, R₁, a₂, R₂, p)

Power mean of the effective resource supplies `X₁ = a₁R₁` and `X₂ = a₂R₂` with
exponent `p.s`: `((X₁^s + X₂^s) / 2)^(1/s)`, or the geometric mean `√(X₁X₂)` for
`s == 0`.
"""
calc_T(a₁, R₁, a₂, R₂, p) = calc_T(a₁ * R₁, a₂ * R₂, p)

function calc_T(X₁, X₂, p)
    (; s) = p
    if s == 0.0
        return sqrt(X₁ * X₂)
    else
        return 2^(-1 / s) * (X₁^s + X₂^s)^(1 / s)
    end
end

"""
    calc_G(X₁, X₂, p)
    calc_G(a₁, R₁, a₂, R₂, p)
    calc_G(χ, R₁, R₂, p)

Saturating response `G = 1 / (α + β / T)`, where `T` is [`calc_T`](@ref).
The `χ` method computes the affinities from the trait and clamps them at zero.
"""
calc_G(X₁, X₂, p) = 1 / (p.α + p.β / calc_T(X₁, X₂, p))
calc_G(a₁, R₁, a₂, R₂, p) = calc_G(a₁ * R₁, a₂ * R₂, p)

function calc_G(χ, R₁, R₂, p)
    a₁, a₂ = calc_a1_a2(χ, p)
    a₁ = max(0.0, a₁)
    a₂ = max(0.0, a₂)
    return calc_G(a₁, R₁, a₂, R₂, p)
end

"""
    calc_C1_C2(X₁, X₂, p)
    calc_C1_C2(a₁, R₁, a₂, R₂, p)

Normalised shares `Cᵢ = Xᵢ^E / (X₁^E + X₂^E)` with `E = exp(s - 1)`.
TODO: describe what `C₁` and `C₂` represent in the model.
"""
calc_C1_C2(a₁, R₁, a₂, R₂, p) = calc_C1_C2(a₁ * R₁, a₂ * R₂, p)

function calc_C1_C2(X₁, X₂, p)
    E = exp(p.s - 1)
    den = X₁^E + X₂^E
    C₁ = X₁^E / den
    C₂ = X₂^E / den
    return C₁, C₂
end

## -----------------------------------------------------------------------------
## Selection gradients
## -----------------------------------------------------------------------------

"""
    calc_local_sel_grad(χ, R₁, R₂, p)

Derivative `∂G/∂χ` at a single grid cell with resource densities `R₁`, `R₂`,
computed by forward-mode automatic differentiation.
"""
function calc_local_sel_grad(χ, R₁, R₂, p)
    return ForwardDiff.derivative(x -> calc_G(x, R₁, R₂, p), χ)
end

"""
    calc_max_sel_grad(y, p)

Largest absolute selection gradient over all species in the state vector `y`.

For each species the gradient is the average of the local gradients
([`calc_local_sel_grad`](@ref)) over the grid, weighted by the squared consumer
density. For a species within `TRAIT_BOUNDARY_TOL` of the trait boundary ±1, only a
gradient pointing back into the interval counts, since the trait cannot evolve
past the boundary.
"""
function calc_max_sel_grad(y, p)
    R₁, R₂, u, χ = get_state_variables_eco_evo(y, p)

    ∂λ∂χ_eff = zero(χ)  # effective gradient of each species
    for i ∈ eachindex(χ, ∂λ∂χ_eff)
        uᵢ = view(u, :, i)
        ∂λ∂χ = sum(uᵢ[j]^2 * calc_local_sel_grad(χ[i], R₁[j], R₂[j], p)
                   for j ∈ eachindex(R₁, R₂, uᵢ))
        ∂λ∂χ /= sum(uᵢⱼ^2 for uᵢⱼ ∈ uᵢ)

        if (1 - abs(χ[i])) >= TRAIT_BOUNDARY_TOL
            ∂λ∂χ_eff[i] = ∂λ∂χ
        elseif χ[i] > 1.0 - TRAIT_BOUNDARY_TOL
            ∂λ∂χ_eff[i] = min(0.0, ∂λ∂χ)
        elseif χ[i] < -1.0 + TRAIT_BOUNDARY_TOL
            ∂λ∂χ_eff[i] = max(0.0, ∂λ∂χ)
        end
    end

    return maximum(abs, ∂λ∂χ_eff)
end

## -----------------------------------------------------------------------------
## Finite-difference variant ("_spar")
## -----------------------------------------------------------------------------
# TODO: confirm whether these are still used. Compared with the functions above, they
# have no `s == 0` branch and do not clamp the affinities at zero.

calc_T_spar(X₁, X₂, p) = 2^(-1 / p.s) * (X₁^p.s + X₂^p.s)^(1 / p.s)
calc_G_spar(X₁, X₂, p) = 1 / (p.α + p.β / calc_T_spar(X₁, X₂, p))
calc_G_spar(a₁, R₁, a₂, R₂, p) = calc_G_spar(a₁ * R₁, a₂ * R₂, p)

function calc_G_spar(χ, R₁, R₂, p)
    a₁, a₂ = calc_a1_a2(χ, p)
    return calc_G_spar(a₁, R₁, a₂, R₂, p)
end

# Central finite-difference estimate of ∂G/∂χ with step `h`.
function calc_local_sel_grad_spar(χ, R₁, R₂, p; h = 1e-6)
    G_a = calc_G_spar(χ - h, R₁, R₂, p)
    G_b = calc_G_spar(χ + h, R₁, R₂, p)
    return (G_b - G_a) / (2h)
end