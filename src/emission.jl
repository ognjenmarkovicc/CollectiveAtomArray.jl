"""
Code for spatial emission.
"""

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
    unit_r = [sin(phi)*sin(theta), cos(phi)*sin(theta),
              cos(theta)]

    return (3/(8*π))*(1-abs2(dot(dipole, unit_r)))
end