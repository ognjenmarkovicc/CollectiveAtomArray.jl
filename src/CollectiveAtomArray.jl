module CollectiveAtomArray

using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra
using Symbolics: substitute, simplify
using QuantumCumulants: QMul

export get_system_operators, get_system_eqs, get_param_substitution

include("system.jl")
include("equations.jl")

end