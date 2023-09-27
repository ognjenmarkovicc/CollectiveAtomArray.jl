"""
Save equations data.
"""

import CollectiveSpins as cs
import CollectiveAtomArray as caa
import CollectiveSpins: SpinCollection

using JLD2
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra

data_folder = "data"

Ns = [2]

for N in Ns
    eqs = caa.get_system_eqs(N);
    jldsave(joinpath(data_folder, "eqs_$(N)_atoms.jld2"); eqs)
end

lattice_spacing = 0.1
dipole = normalize([0.0, 0.0, 1.0])
time_start = 0.0
time_end = 5.0

for N in Ns
    system = cs.SpinCollection(cs.geometry.chain(lattice_spacing, N), 
    dipole, 1.); # Γ₀=1.

    eqsfile = load(joinpath(data_folder, "eqs_$(N)_atoms.jld2"));

    odesys, u0, param_subs = caa.get_system_ode_from_eqs(eqsfile["eqs"], system);

    prob = ODEProblem(odesys,u0,(time_start, time_end),param_subs);

    sol = solve(prob,RK4());

end