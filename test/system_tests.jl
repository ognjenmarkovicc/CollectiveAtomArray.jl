"""
Test collective system setup.
"""

@testset "interface with CollectiveSpins" begin

    N=2

    system = SpinCollection(geometry.chain(0.1, N), 
                            [0., 0., 1.], 1.); 

    pos = get_system_pos(system)

    @test length(system) == N
    @test size(pos) == (3, N)
    @test pos[:, 1] == system.spins[1].position
end