module CollectiveAtomArray

using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra
using Symbolics: substitute, simplify
using QuantumCumulants: QMul
import CollectiveSpins: SpinCollection, interaction

export get_system_sigma, get_system_eqs, 
get_param_substitution, get_system_ode, get_emission_factor,
direction_vec

include("system.jl")
include("solution.jl")
include("emission.jl")

end