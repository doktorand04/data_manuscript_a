"""
    ADRD

TODO: one-paragraph description of the model (consumers with a trait trade-off
competing for two diffusing resources on a periodic grid) and what the module
computes (eco-evolutionary trajectories to an evolutionarily stable state).

The main entry points are [`biodiversity_generator_2`](@ref) (start from a single
species) and [`biodiversity_generator_from_mid`](@ref) (resume from a saved state).
"""
module ADRD

# Linear algebra, automatic differentiation, updating immutable structs
using LinearAlgebra
using ForwardDiff
using Accessors: @set

# ODE solving
using OrdinaryDiffEq                        # ODEProblem, solve, KenCarp4, callbacks
using DiffEqCallbacks: TerminateSteadyState
using Sundials: CVODE_BDF
using LinearSolve: KrylovJL_GMRES
import SciMLBase                            # SciMLBase.AutoSpecialize

# Preallocated work arrays that are compatible with automatic differentiation
import PreallocationTools

# Reading and writing .xlsx files
import DataFrames
import XLSX

include("spectral.jl")    # spectral grid, differentiation matrices, interpolation
include("parameters.jl")  # EcoEvoParameters and ODECache
include("fitness.jl")     # affinities, growth function, selection gradients
include("dynamics.jl")    # ODE right-hand sides, solver wrappers, state vectors
include("invasion.jl")    # invasion fitness and invader selection
include("extinction.jl")  # removal of extinct species
include("simulation.jl")  # eco-evolutionary driver and result I/O

end # module
