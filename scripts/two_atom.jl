import CollectiveSpins as cs
import CollectiveAtomArray as caa

using Plots
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra

# parameters
N = 2 # number of atoms
lattice_spacing = 0.1
dipole = [0., 0., 1.]

# system operators
h, σ = caa.get_system_operators(N)

eqs = caa.get_system_eqs(N, order=2)

# Generate the ODESystem
@named sys = ODESystem(eqs)

# Generate coupling matrices
system = cs.SpinCollection(cs.geometry.chain(lattice_spacing, N), 
                           dipole, 1.) # Γ₀=1.

Γmat = cs.interaction.GammaMatrix(system)
Jmat = cs.interaction.OmegaMatrix(system)

# variable replacements
p = caa.get_param_substitution(Jmat, Γmat)

# Create ODEProblem

# initial state -- all atoms in the excited states
u0 = convert(Array{ComplexF64}, map(caa.get_initial_all_excited,
                                    eqs.operators))
prob = ODEProblem(sys,u0,(0.0,5.0),p)

# Solve
sol = solve(prob,RK4())

# Plot
graph = plot(sol.t, real.(sol[σ(:e,:e,1)]), label="Excited state in atom 1",
             xlabel="γt", ylabel="Excited state population")
plot!(graph, sol.t, real.(sol[σ(:e,:e,N)]), label="End of chain", leg=1)


# list of the decay operators
decay_rate = sum([Γmat[i, j]*real(sol[σ(:e, :g, i)*σ(:g, :e, j)]) for i=1:N for j=1:N])

# Plot decay rate
graph = plot(sol.t, decay_rate, xlabel="γt", ylabel="Total decay rate")