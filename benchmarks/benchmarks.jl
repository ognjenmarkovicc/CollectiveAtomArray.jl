import CollectiveSpins as cs
import CollectiveAtomArray as caa
using ModelingToolkit, OrdinaryDiffEq

using BenchmarkTools
using Suppressor

function benchmark_solution()
    # parameters
    N = 4 # number of atoms
    lattice_spacing = 0.1
    dipole = normalize([0., 0., 1.])
    time_start = 0.0
    time_end = 5.0 

    system = cs.SpinCollection(cs.geometry.chain(lattice_spacing, N), 
                            dipole, 1.);

    # benchmark getting problem
    println("Benchmarking CollectiveAtomArray.get_system_ode")
    ben = @benchmark odesys, u0, param_subs = caa.get_system_ode($system)
    show(stdout, "text/plain", ben)

    # benchmark getting equations
    println("Benchmarking CollectiveAtomArray.get_system_eqs")
    ben = @benchmark eqs = caa.get_system_eqs($N, order=2)
    # print the benchmark (for some reason only the last one is
    #printed in vscode)
    show(stdout, "text/plain", ben)

    odesys, u0, param_subs = caa.get_system_ode(system)
    prob = ODEProblem(odesys,u0,(time_start,time_end),param_subs);
    # benchmark solving the equation solve
    println("Benchmarking diff equation solve")
    # suppress retcode warning
    @suppress begin
        @benchmark sol = solve($prob,RK4())
    end
end

benchmark_solution()