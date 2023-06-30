import CollectiveSpins as cs
import CollectiveAtomArray as caa

using Plots
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra

# parameters
N = 10 # number of atoms
lattice_spacing = 0.1
dipole = [0., 0., 1.]
time_start = 0.0
time_end = 5.0 

system = cs.SpinCollection(cs.geometry.chain(lattice_spacing, N), 
                           dipole, 1.) # Γ₀=1.

odesys, u0, param_subs = caa.get_system_ode(system)

prob = ODEProblem(odesys,u0,(time_start,time_end),param_subs)

# Solve
sol = solve(prob,RK4())

# Plot
_, σ = caa.get_system_operators(N)
Γmat = cs.interaction.GammaMatrix(system)

# plot excitation fraction
exc_frac = sum([real(sol[σ(:e, :e, i)]) for i=1:N]) 

g1 = plot(sol.t, exc_frac/N, xlabel="γt", ylabel="Excitation fraction")

# list of the decay operators
decay_rate = sum([Γmat[i, j]*real(sol[σ(:e, :g, i)*σ(:g, :e, j)]) for i=1:N for j=1:N])

# Plot decay rate
g2 = plot(sol.t, decay_rate/N, 
          xlabel="γt", ylabel="Decay rate/(N Γ₀)")

g = plot(g1, g2, xscale=:log10, xlim=(1e-2, 5),
         layout=(2, 1), legend=false)