"""
Example script to 
replicate Figure 1b, PRR 4, 023207 (2022)
"""

import CollectiveSpins as cs
import Statistics as st
using Plots

dipole = [0., 0., 1.]
N = 25 # atom number
Γ₀ = 1. 

lattice_spacings = collect(0.1:0.01:2)
Γ_var = zeros(length(lattice_spacings))
for (idx, lattice_spacing) in enumerate(lattice_spacings)
    system = cs.SpinCollection(cs.geometry.chain(lattice_spacing, N), 
                               dipole, Γ₀)

    Γmat = cs.interaction.GammaMatrix(system)

    eigen_vals = eigen(Γmat).values
    Γ_var[idx] = st.var(eigen_vals)
end

graph = plot(lattice_spacings, Γ_var, xlabel="Lattice spacing",
             ylabel="var(Γ)", ylim=(0, 2),
             xlim=(0, 2))
