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

"""
get_mean_intensity(θ::Float64, ϕ::Float64,
    mean_correlation::Union{Matrix{Float64}, Matrix{Vector{Float64}}},,
    system::SpinCollection)

Get the mean intensity radiated in direction (θ,ϕ). mean_correlation can be
a matrix of numbers of shape (N,N) or a matrix of vectors, where each vector
represents correlations for different time points.

mean_correlation is returned by get_mean_correlation
"""
function get_mean_intensity(θ::Float64, ϕ::Float64,
    mean_correlation::Union{Matrix{Float64}, Matrix{Vector{Float64}}},
    system::SpinCollection)
    # pos: size (3, N) matrix of spin positions
    pos =  get_system_pos(system)

    # all exp(-kr̂⋅(rᵢ-rⱼ)) factors in a (N, N) matrix
    compl_factors = get_compl_factors(θ, ϕ, pos)

    # the dipole emission factor
    # assuming all the atoms have the same polarization
    dipole = system.polarizations[1]
    dip_factor = get_emission_factor(θ, ϕ, dipole)

    return dip_factor*sum(compl_factors .* mean_correlation)
end
