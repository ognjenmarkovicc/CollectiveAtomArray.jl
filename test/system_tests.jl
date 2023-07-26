"""
Test collective system setup.
"""

@testset "interface with CollectiveSpins" begin

    N=2

    system = SpinCollection(geometry.chain(0.1, N), 
                            [0., 0., 1.], 1.); 

    @test length(system) == N
end