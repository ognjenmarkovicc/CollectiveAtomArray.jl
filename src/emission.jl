"""
Code for spatial emission.
"""

"""
    direction_vec(theta::Real, phi::Real)

Get a unit vector in direction (θ, ϕ)

# Arguments
* `theta::Real`: Polar angle.
* `phi::Real`: Azimuthal angle.
"""
function direction_vec(theta::Real, phi::Real)
    return [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
end

"""
    get_emission_factor(theta::Real, phi::Real, dipole::Vector)

Get the dipole emission factor. 

# Arguments
* `theta::Real`: Polar angle.
* `phi::Real`: Azimuthal angle.
* `dipole::Vector`: Dipole vector
"""
function get_emission_factor(theta::Real, phi::Real, dipole::Vector)
    dipole = normalize(dipole) # make dipole a unit vector
    unit_r = direction_vec(theta, phi)

    return (3/(8*π))*(1-abs2(dot(dipole, unit_r)))
end

"""
get_compl_factors(theta::Real, phi::Real, pos::Matrix{Real})

Get the complex exp(-kr̂⋅(rᵢ-rⱼ)) factors in a matrix,
where r̂ is the unit direction vector and rᵢ, rⱼ are
the spin position vectors.
"""
function get_compl_factors(theta::Real, phi::Real, pos::Matrix{Float64})
    unit_r = direction_vec(theta, phi)
    # @assert size(pos)[1] == 3, "First dimension of post must have length 3"
    (_, N) = size(pos)
    return [exp(-2*im*π*dot(unit_r, pos[:, i] - pos[:, j])) for i=1:N, j=1:N]
end
