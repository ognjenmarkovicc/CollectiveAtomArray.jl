module CollectiveAtomArray

using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra
using Symbolics: substitute, simplify
using QuantumCumulants: QMul
import CollectiveSpins: SpinCollection, interaction

export get_system_operators, get_system_eqs, get_param_substitution, get_system_ode

include("system.jl")
include("equations.jl")

end