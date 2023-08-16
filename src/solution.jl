"""
Code for setting up the system equations and solving them.
"""

"""Get dissipative coupling parameter"""
Γparam(i, j) = cnumber("Γ_{$i, $j}")
"""Get coherent coupling parameter"""
Jparam(i, j) = cnumber("J_{$i, $j}")

"""
    get_system_eqs(N::Integer, order::Integer=2)

# Arguments
* `N::Integer`: Number of atoms in the system.
* `order::Integer`: meanfield equation order
"""
function get_system_eqs(N::Integer; order::Integer=2)
    # define dissipative and coherent coupling
    # as a vector of vectors of cnumbers
    Γ = [[Γparam(i, j) for j in 1:N]; for i in 1:N]

    J = [[Jparam(i, j) for j in 1:N]; for i in 1:N]

    # Symbol(:atom, i) results in symbol "atom$(i)$,
    # ie if i=1, then it is atom1
    σ = get_system_sigma(N)

    # coherent Hamiltonian  
    H = sum(J[i][j]*σ(:e, :g, i)*σ(:g, :e, j) for i in 1:N for j in 1:N)

    # jump operators
    # build the list of correlated jump operators
    corr_c_ops = [sum(Γ[i][j]*σ(:g,:e,j) for j in 1:N) for i=1:N]

    # hermitian conjugates are just the uncorrelated operators
    c_ops_dagg = [σ(:e,:g,j) for j in 1:N]

    relevant_ops = get_nonzero_operators(N, 
                                         order=order)

    # depending on c_ops and zero and non-zero terms
    # we will have different number of operators
    eqs = meanfield(relevant_ops,H,corr_c_ops;
                    Jdagger=c_ops_dagg,
                    rates=[1. for i=1:N],
                    order=order,
                    simplify=false)
    # find states which are on the equation rhs
    missing_states = find_missing(eqs)
    # set them to 0
    missing_subs = Dict(missing_states .=> 0)
    eqs = substitute(eqs, missing_subs)

    return eqs
end

"""
    get_param_substitution(Jmat::Matrix,
                           Γmat::Matrix)

Replace symbolic coupling parameters with actual coupling
values taken from a provided matrices of values.
"""
function get_param_substitution(Jmat::Array{<:AbstractFloat, 2},
                                Γmat::Array{<:AbstractFloat, 2})
    Jsize = size(Jmat)
    Γsize = size(Γmat)

    Jmat = reduce(vcat, Jmat)
    Γmat = reduce(vcat, Γmat)

    # flattened coupling parameter vectors
    Jflat = [Jparam(i, j) for i=1:Jsize[1] for j=1:Jsize[2]]
    Γflat = [Γparam(i, j) for i=1:Γsize[1] for j=1:Γsize[2]]
    return [Γflat .=> Γmat; Jflat.=>Jmat;]
end

"""
    get_system_ode(system <: SpinCollection, order::Integer=2)

Get the system odesys (<: ModelingToolkit.AbstractODESystem),
initial conditions and system parameter substitutions (i.e. Γᵢⱼ -> specific value))
"""
function get_system_ode(system::SpinCollection,
                        order::Integer=2)

    N = length(system.spins)

    eqs = get_system_eqs(N, order=order)

    # Generate the ODESystem
    @named odesys = ODESystem(eqs)

    # ensure Γmat is real, as it is from definition
    # e.g. PRR 5, 013091 (2023), Eq. (4)
    Γmat = real(interaction.GammaMatrix(system))
    Jmat = interaction.OmegaMatrix(system)

    # variable replacements
    param_subs = get_param_substitution(Jmat, Γmat)

    # Create ODEProblem
    # initial state -- all atoms in the excited states
    u0 = convert(Array{ComplexF64}, map(get_initial_all_excited,
                                        eqs.operators))

    return odesys, u0, param_subs
end

"""
    get_mean_excitation(system::SpinCollection, sol::OrdinaryDiffEq.ODESolution)

    Get <σ_iee> for each i ∈ 1:N. Returns a length N vector of vectors,
    each of length(sol.t).
"""
function get_mean_excitation(system::SpinCollection, sol::OrdinaryDiffEq.ODESolution)
    σ = get_system_sigma(system)
    return [real(sol[σ(:e, :e, i)]) for i=1:length(system)]
end

"""
    get_mean_correlation(system::SpinCollection, sol::OrdinaryDiffEq.ODESolution)

    Get <σ_ieg σ_jge> for each i, j ∈ 1:N. Returns a (N, N) matrix of vectors,
    each of length(sol.t).
"""
function get_mean_correlation(system::SpinCollection, sol::OrdinaryDiffEq.ODESolution)
    σ = get_system_sigma(system)
    N = length(system)
    return [real(sol[σ(:e, :g, i)*σ(:g, :e, j)]) for i=1:N, j=1:N]
end

"""
    get_decay_rate(system::SpinCollection, sol::OrdinaryDiffEq.ODESolution)

    Get decay rate from the solution to the excitation dynamics.
"""
function get_decay_rate(system::SpinCollection, sol::OrdinaryDiffEq.ODESolution)
    Γmat = real(interaction.GammaMatrix(system))
    mean_correlation = get_mean_correlation(system, sol);

    # decay rate is the sum of <σ_ieg σ_jge>
    # correlations weighted by the Γ matrix
    # PRR 5, 013091 (2023), Eq. (8)
    return sum(Γmat.*mean_correlation)
end