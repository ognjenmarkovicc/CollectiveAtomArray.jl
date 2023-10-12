using CollectiveAtomArray
using LinearAlgebra
using Test
using CollectiveSpins
using Cubature

@testset "CollectiveAtomArray.jl" begin
    @testset "Equation tests" begin
        include("equation_tests.jl")
    end
    @testset "Spatial emission tests" begin
        include("emission_tests.jl")
    end
    @testset "System tests" begin
        include("system_tests.jl")
    end
end
