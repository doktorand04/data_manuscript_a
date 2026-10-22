# Invasion fitness and invader selection.
#
# This file is meant to be `include`d from ADRD.jl, which should load LinearAlgebra
# (for `diag` and `eigvals`).

# Candidate invaders closer than this to a resident trait are ignored.
const MIN_INVADER_TRAIT_DISTANCE = 1e-2
# A candidate is accepted as an invader only if its invasion fitness is at least this.
const INVASION_FITNESS_TOL = 1e-8

"""
    calc_inv_fit!(F, δΔ_diag, χₘ, R₁ₛ, R₂ₛ, p)

Invasion fitness of a rare invader with trait `χₘ` in the resource environment
`R₁ₛ`, `R₂ₛ`: the largest eigenvalue of `δΔ + diag(G - μ)`, with `μ = 1/α - 1`
(the linearised consumer dynamics).

`F` must be a copy of `p.δΔ` and `δΔ_diag` a copy of its diagonal. Only the diagonal
of `F` is overwritten, so the same `F` can be reused for many traits.
"""
function calc_inv_fit!(F, δΔ_diag, χₘ, R₁ₛ, R₂ₛ, p)
    for j ∈ eachindex(R₁ₛ, R₂ₛ, δΔ_diag)
        F[j, j] = δΔ_diag[j] + calc_G(χₘ, R₁ₛ[j], R₂ₛ[j], p) + 1 - 1 / p.α
    end
    return maximum(real, eigvals(F))
end

"""
    calc_inv_fit(χₘ, R₁ₛ, R₂ₛ, p)

Invasion fitness of a single trait `χₘ`; see [`calc_inv_fit!`](@ref).
"""
function calc_inv_fit(χₘ, R₁ₛ, R₂ₛ, p)
    F = copy(p.δΔ)
    return calc_inv_fit!(F, diag(F), χₘ, R₁ₛ, R₂ₛ, p)
end

"""
    calc_fitness_landscape(χₘ, R₁ₛ, R₂ₛ, p)
    calc_fitness_landscape(y, p)

Invasion fitness for each trait in `χₘ`, given resource densities `R₁ₛ`, `R₂ₛ`.
The second method evaluates the landscape at the resident traits of the
eco-evolutionary state vector `y`.
"""
function calc_fitness_landscape(χₘ, R₁ₛ, R₂ₛ, p)
    F = copy(p.δΔ)
    δΔ_diag = diag(F)  # a copy: the diagonal of F is overwritten below
    λds = zeros(length(χₘ))
    for i ∈ eachindex(χₘ)
        λds[i] = calc_inv_fit!(F, δΔ_diag, χₘ[i], R₁ₛ, R₂ₛ, p)
    end
    return λds
end

function calc_fitness_landscape(y, p)
    R₁ₛ, R₂ₛ, _, χₘ = get_state_variables_eco_evo(y, p)
    return calc_fitness_landscape(χₘ, R₁ₛ, R₂ₛ, p)
end

"""
    calc_self_inv_fits(y_eco, χ_ress, p)

Invasion fitness of each trait in `χ_ress` in the resource environment of the
ecological state vector `y_eco`.
"""
function calc_self_inv_fits(y_eco, χ_ress, p)
    R₁ₛ, R₂ₛ, _ = get_state_variables_eco(y_eco, p)
    return calc_fitness_landscape(χ_ress, R₁ₛ, R₂ₛ, p)
end

"""
    calc_maximally_fit_invader_trait(χₘ, y, p) -> Vector{Float64}

Pick the best invader from the candidate traits `χₘ`, given the eco-evolutionary
state `y`.

Candidates within `MIN_INVADER_TRAIT_DISTANCE` of a resident trait are excluded. Of
the rest, the one with the highest invasion fitness is returned as a one-element
vector, provided its fitness is at least `INVASION_FITNESS_TOL`. Otherwise an empty
vector is returned: no invader can establish, so the community is evolutionarily
stable.
"""
function calc_maximally_fit_invader_trait(χₘ, y, p)
    R₁ₛ, R₂ₛ, _, χₛ = get_state_variables_eco_evo(y, p)

    # Restrict the candidate pool to traits that are not close to a resident
    mindiff = vec(minimum(abs.(χₛ' .- χₘ); dims = 2))
    candidates = χₘ[findall(>(MIN_INVADER_TRAIT_DISTANCE), mindiff)]
    isempty(candidates) && return Float64[]

    λds = calc_fitness_landscape(candidates, R₁ₛ, R₂ₛ, p)
    λ_best, i_best = findmax(λds)
    λ_best < INVASION_FITNESS_TOL && return Float64[]

    return [candidates[i_best]]
end