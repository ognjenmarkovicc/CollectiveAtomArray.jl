using CollectiveAtomArray
using Test

@testset "CollectiveAtomArray.jl" begin
    @testset "Equation tests" begin
        include("equation_tests.jl")
    end
end
