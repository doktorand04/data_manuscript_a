"""
    ADRD

The `ADRD` module implements the eco-evolutionary model of two resources and 
multiple consumer species, as described in the manuscript. It provides functions
for simulating the dynamics of the system, computing invasion fitness, and 
generating biodiversity through adaptive dynamics.

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
