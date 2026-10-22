# Fourier spectral grid, differentiation matrices and interpolation on a periodic
# domain.
#
# This file is meant to be `include`d from ADRD.jl, which should load
# LinearAlgebra (for `Diagonal`).

function _check_even(N)
    (N > 0 && iseven(N)) ||
        throw(ArgumentError("Only even N permissible (got N = $N)"))
    return nothing
end

"""
    calc_spectral_x(N; L = 2π)

Return the `N` equispaced points `h, 2h, …, L` (with `h = L/N`) of a periodic
domain of length `L`. `N` must be even.
"""
function calc_spectral_x(N; L = 2π)
    _check_even(N)
    h = L / N
    return h * (1:N)
end

"""
    calc_spectral_x_and_D₂(N; L = 2π) -> (x, D₂)

Return the grid `x` from [`calc_spectral_x`](@ref) and the Fourier spectral
second-derivative matrix `D₂` for periodic functions sampled on it. `N` must be even.

This is the standard periodic differentiation matrix (e.g. Trefethen,
*Spectral Methods in MATLAB*, SIAM 2000, ch. 3), rescaled from period `2π` to `L`.
"""
function calc_spectral_x_and_D₂(N; L = 2π)
    x = calc_spectral_x(N; L)

    # D₂ is symmetric and circulant, so entry (i, j) depends only on |i - j|.
    column = [k == 0 ? -N^2 / 12 - 1 / 6 : -0.5 * (-1)^k * csc((π / N) * k)^2
              for k in 0:N-1]
    D₂ = [column[abs(i - j) + 1] for i in 1:N, j in 1:N]
    D₂ .*= (2π / L)^2

    return x, D₂
end

"""
    calc_spectral_D22(N; L = 2π)
    calc_spectral_D22(D₂::AbstractMatrix)

Two-dimensional spectral Laplacian on an `N × N` periodic grid, as a dense
`N² × N²` matrix. Built from the 1D second-derivative matrix `D₂` as
`I ⊗ D₂ + D₂ ⊗ I`.
"""
function calc_spectral_D22(N; L = 2π)
    _, D₂ = calc_spectral_x_and_D₂(N; L)
    return calc_spectral_D22(D₂)
end

function calc_spectral_D22(D₂::AbstractMatrix)
    id  = Diagonal(ones(size(D₂, 1)))
    D₂₂ = Matrix(kron(id, D₂) + kron(D₂, id))

    # Reset each diagonal entry so that every row sums to exactly zero
    # (a constant field then has exactly zero Laplacian).
    for i in axes(D₂₂, 1)
        D₂₂[i, i] = 0
        D₂₂[i, i] = -sum(view(D₂₂, i, :))
    end

    return D₂₂
end

"""
    spectral_interpolate(A, N_new::Integer; L = 2π)
    spectral_interpolate(A, xs_new; L = 2π)

Fourier-interpolate the periodic field `A` (an `N × N` array sampled on the grid of
[`calc_spectral_x`](@ref), `N` even) onto a new `m × m` grid, using the periodic sinc
kernel. The new grid is given either as a number of points `N_new` (points
`L/N_new, 2L/N_new, …, L`) or as an explicit vector of coordinates `xs_new`.
"""
function spectral_interpolate(A, N_new::Integer; L = 2π)
    xs_new = (1/N_new:1/N_new:1) * L
    return spectral_interpolate(A, xs_new; L)
end

function spectral_interpolate(A, xs_new; L = 2π)
    N = size(A, 1)
    size(A, 2) == N ||
        throw(DimensionMismatch("A must be square, got size $(size(A))"))
    _check_even(N)

    x   = calc_spectral_x(N; L)
    Lsf = 2π / L

    Xd = xs_new .- x'   # Xd[i, j] = xs_new[i] - x[j]
    S  = @. sin(Lsf * Xd * N * 0.5) / (N * tan(Lsf * Xd * 0.5))
    S[iszero.(Xd)] .= 1  # limit of the kernel at zero offset

    return S * A * S'
end