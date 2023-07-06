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
    return [sin(phi)*sin(theta), cos(phi)*sin(theta), cos(theta)]
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