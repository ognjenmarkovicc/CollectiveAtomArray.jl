@testset "System operators" begin

    N=2

    h, σ = get_system_operators(N)

    eqs = get_system_eqs(N, order=2)

    operators = Set(eqs.operators)

    @test σ(:e, :e, 1) in operators
    @test σ(:e, :e, 2) in operators
    @test σ(:e, :g, 1)*σ(:g, :e, 2) in operators
    @test σ(:e, :e, 1)*σ(:e, :e, 2) in operators

end