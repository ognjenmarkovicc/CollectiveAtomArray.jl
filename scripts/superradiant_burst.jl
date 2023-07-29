"""
Example script to 
replicate Figure 2 (a) and (c), PRR 5, 013091 (2023)
Currently only replicates the second order cumulant
expansion simulation.
"""

import CollectiveSpins as cs
import CollectiveAtomArray as caa

using Plots
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra

# parameters
N = 10 # number of atoms
lattice_spacing = 0.1
dipole = normalize([0., 0., 1.])
time_start = 0.0
time_end = 5.0 

system = cs.SpinCollection(cs.geometry.chain(lattice_spacing, N), 
                           dipole, 1.); # Γ₀=1.

odesys, u0, param_subs = caa.get_system_ode(system)

prob = ODEProblem(odesys,u0,(time_start,time_end),param_subs)

# Solve
sol = solve(prob,RK4())

# Plot
# plot excitation fraction
exc_frac = sum(get_mean_excitation(system, sol)); 

g1 = plot(sol.t, exc_frac/N, xlabel="γt", ylabel="Excitation fraction");

decay_rate = get_decay_rate(system, sol);

# Plot decay rate
g2 = plot(sol.t, decay_rate/N, 
          xlabel="γt", ylabel="Decay rate/(N Γ₀)")

g = plot(g1, g2, xscale=:log10, xlim=(1e-2, 5),
         layout=(2, 1), legend=false)