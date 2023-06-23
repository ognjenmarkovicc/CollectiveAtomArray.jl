@testset "System operators N=2" begin

    N=2

    h, σ = get_system_operators(N)

    eqs = get_system_eqs(N, order=2)

    operators = Set(eqs.operators)

    @test σ(:e, :e, 1) in operators
    @test σ(:e, :e, 2) in operators
    @test σ(:e, :g, 1)*σ(:g, :e, 2) in operators
    @test σ(:e, :e, 1)*σ(:e, :e, 2) in operators

end

@testset "System operators N=3" begin

    N=3

    h, σ = get_system_operators(N)

    eqs = get_system_eqs(N, order=2)

    operators = Set(eqs.operators)
    @show operators
    @test σ(:e, :e, 1) in operators
    @test σ(:e, :e, 2) in operators
    @test σ(:e, :g, 1)*σ(:g, :e, 2) in operators
    @test σ(:e, :e, 1)*σ(:e, :e, 2) in operators
    @test σ(:g, :e, 1) ∉ operators
    @test σ(:e, :g, 1) ∉ operators


end