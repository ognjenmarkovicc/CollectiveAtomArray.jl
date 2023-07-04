@testset "spatial emission tests" begin
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