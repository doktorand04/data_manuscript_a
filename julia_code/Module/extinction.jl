# Removal of extinct species from the ecological state.
#
# This file is meant to be `include`d from ADRD.jl, which should load a package that
# provides `@set` (Accessors or Setfield).

# A species whose maximum density over the grid is at or below this is considered
# extinct. (Not to be confused with `ODE_EXTINCTION_LIMIT` in dynamics.jl, the density
# floor that keeps the ODE solution non-negative.)
const EXTINCTION_DENSITY_THRESHOLD = 1e-6

"""
    extinction(χ, y, p) -> (χ, y, p)

Remove the extinct species from the ecological state vector `y`: those whose maximum
density over the grid is at most `EXTINCTION_DENSITY_THRESHOLD`.

Returns the traits of the surviving species, the new state vector, and a copy of `p`
with `S` updated. The input `y` is not modified.
"""
function extinction(χ, y, p)
    R₁, R₂, U = get_state_variables_eco(y, p)

    extinct = vec(maximum(U; dims = 1)) .<= EXTINCTION_DENSITY_THRESHOLD
    keep = findall(.!extinct)
    isempty(keep) &&
        error("All species went extinct (maximum density ≤ $EXTINCTION_DENSITY_THRESHOLD)")

    χ = χ[keep]
    U = U[:, keep]
    p = @set p.S = length(keep)
    y = get_state_vector_eco(R₁, R₂, U)

    return χ, y, p
end