@testset "dipole emission factor" begin
    """
    Test the numerical value of the dipole emission factor.
    E.g. for linear polarization, it should be 
    (3/(8π))*sin²(theta)
    """
    dipole = normalize([0, 0, 1])

    @test get_emission_factor(0, 0, dipole) ≈ 0
    @test get_emission_factor(π/2, 0, dipole) ≈ 3/(8*π)

    dipole = normalize([1, 1im, 0])

    @test get_emission_factor(0, 0, dipole) ≈ 3/(8*π)
    @test get_emission_factor(π/2, 0, dipole) ≈ 3/(16*π)
end

@testset "unit direction vec" begin
    @test isapprox(direction_vec(0., 0.), [0., 0., 1.], atol=1e-11)
    @test isapprox(direction_vec(π/2, 0.), [1., 0., 0.], atol=1e-11)
end

@testset "intensity calculation" begin
    test_pos = [[0.1, 0, 0] [0.2, 0, 0]]
    test_compl_factor = Matrix{ComplexF64}([[1, 1] [1, 1]])
    @test get_compl_factors(0., 0., test_pos) ≈ test_compl_factor

    # correlations in an uncorrelated system
    N=4
    mean_correlation_test = diagm(ones(4))
    dipole = normalize([0, 0, 1])
    system = SpinCollection(geometry.chain(1.0, 4), dipole, 1.)
    test_emission_factor = get_emission_factor(0.,0.,dipole)
    @test get_mean_intensity(0.,0.,mean_correlation_test, system)/N ≈ test_emission_factor

    # test the integral over space of mean intensity
    function integrand(vars)
        θ, ϕ = vars
        return sin(θ)*real(get_mean_intensity(θ, ϕ, mean_correlation_test,system))
    end
    result, _ = hcubature(integrand, [0, 0], [π, 2π])
    @test result ≈ N
end