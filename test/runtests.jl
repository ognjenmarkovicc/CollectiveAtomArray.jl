using CollectiveAtomArray
using LinearAlgebra
using Test

@testset "CollectiveAtomArray.jl" begin
    @testset "Equation tests" begin
        include("equation_tests.jl")
    end
    @testset "Spatial emission tests" begin
        include("emission_tests.jl")
    end
end
