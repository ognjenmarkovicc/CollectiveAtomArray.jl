"""
Example script to integrate the total emission rate over observation different angles.
""" 

import CollectiveSpins as cs
import CollectiveAtomArray as caa

using Plots
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra
using Printf
using Cubature

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

# do integration of a vector function
# as get_mean_intensity returns a vector when
# the correlation is evaluated at multiple points
function integrand(vars, v)
    θ, ϕ = vars
    v[:] = sin(θ)*real(caa.get_mean_intensity(θ, ϕ, mean_correlation,system))
end

# decay rate only in a range of NAs
NAmin=0.35
NAmax=0.9

θmin=asin(NAmin)
θmax=asin(NAmax)

integrated_intensity, _ = hcubature(length(sol.t), integrand,
                                    [θmin, 0], [θmax, 2π],reltol=1e-6)

# plot the decay rate vs. time
p = plot(sol.t, integrated_intensity/N, 
          xlabel="γt", ylabel="Decay rate/(N Γ₀)",
          xlim=(1e-2, 1),
          xscale=:log10,
          label="N=$(N) integrated",
          title="Integrated from NAmin=$(NAmin) to NAmax=$(NAmax)")
