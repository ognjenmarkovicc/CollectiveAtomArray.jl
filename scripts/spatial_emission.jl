"""
Example script to 
replicate Figure 3 (a) top left of, PRR 125, 263601 (2020),
which represents the spatial distribution of emitted 
photons from an atomic array for different times and emission angles.
""" 

import CollectiveSpins as cs
import CollectiveAtomArray as caa

using Plots
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra
using Printf

# parameters
N = 10 # number of atoms
lattice_spacing = 0.1
dipole = normalize([0., 1, 1im])
time_start = 0.0
time_end = 20.0 

system = cs.SpinCollection(cs.geometry.chain(lattice_spacing, N), 
                           dipole, 1.); # Γ₀=1.

odesys, u0, param_subs = caa.get_system_ode(system)

prob = ODEProblem(odesys,u0,(time_start,time_end),param_subs)

# Solve
sol = solve(prob,RK4())

# Plot spatial distribution
mean_correlation = caa.get_mean_correlation(system, sol);

# pos: shape (3, N) matrix of spin positions
pos =  caa.get_system_pos(system)

ϕ = 0.

thetas = collect(-π/2:0.1:π/2);
mean_intensity = zeros(ComplexF64, length(thetas), length(sol.t));

for (idx, θ) in enumerate(thetas)

    mean_intensity[idx, :] = caa.get_mean_intensity(θ,ϕ,mean_correlation,system)

end

# plot only the selected times
selected_times = collect(0:2:20)
time_indices = [argmin(abs.(sol.t .- t)) for t in selected_times]

# normalize by the mean intensity at time 0
norm_intensity = real(mean_intensity)[:, time_indices]/real(mean_intensity[1, 1])

labels = [@sprintf "%.1f Γt" t 
          for t in sol.t[time_indices]]

g = plot(thetas, norm_intensity, 
          xlabel="θ", ylabel="Mean normalized intensity", legend=true,
          yscale=:log10, label=reshape(labels, (1, :)))        