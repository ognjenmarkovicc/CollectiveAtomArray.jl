module CollectiveAtomArray

using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using LinearAlgebra
using Symbolics: substitute, simplify
using QuantumCumulants: QMul

export get_system_operators, get_system_eqs, get_param_substitution

"""
    get_system_operators(N)
    
Get the Hilbert space and the sigma operator
for a system of N two-level atoms.
"""
function get_system_operators(N)
    # Hilbert space of the N two-level atoms
    # in this case NLevelSpace is a two-level space
    h = ⊗([NLevelSpace(Symbol(:atom, i), (:g, :e)) for i=1:N]...)

    # operator for the transition from level j to level i for atom k
    # defined as a function
    σ(i,j,k) = Transition(h,Symbol(:σ_,k),i,j,k)

    return h, σ
end

op_is_ee(op::Transition) = op.i == :e && op.j == :e

"""
    get_initial_all_excited(op::Transition)

Get the initial value of the <op> when all atoms are excited.
"""
function get_initial_all_excited(op::Transition)
    # transition between :e and :e
    if op_is_ee(op)
        return 1.
    else
        throw(ArgumentError("$op transition not \\sigma_{ee}"))
    end
end

"""
    get_initial_all_excited(ops::QMul)

Get the initial value the multiplication of operators,
when all atoms are excited.

# Examples

    ⟨σ_1ee σ_2ee⟩=1
    ⟨σ_1ge σ_2eg⟩=0
"""
function get_initial_all_excited(ops::QMul)

    # assume only get two operators
    # if one of the operators is not ee then
    if length(ops.args_nc) == 2
        op1 = ops.args_nc[1]
        op2 = ops.args_nc[2]

        # only if both operators are ee
        return op_is_ee(op1) && op_is_ee(op2) ? 1. : 0.
    else
        throw(ArgumentError("Number of operators in ops needs to be 2"))
    end
end

"""Get dissipative coupling parameter"""
Γparam(i, j) = cnumber("Γ_{$i, $j}")
"""Get coherent coupling parameter"""
Jparam(i, j) = cnumber("J_{$i, $j}")


"""
    get_nonzero_operators(N::Integer, order::Integer)

Get operators which are non-zero for the specific
order of the meanfield equations. This assumes
that the initial system state is incoherent, with 
an arbitrary number of excited atoms.

# Arguments
* `N::Integer`: Number of atoms in the system.
* `order::Integer`: meanfield equation order
"""
function get_nonzero_operators(N::Integer; order::Integer)

    _, σ = get_system_operators(N)

    if order == 2
        relevant_ops = [σ(:e, :g, i)*σ(:g, :e, j) for i in 1:N for j in 1:N if j>=i]
        append!(relevant_ops, [σ(:e, :e, i)*σ(:e, :e, j) for i in 1:N for j in 1:N if j>i])

        return relevant_ops 
    else
        throw(ArgumentError("order in get_nonzero_operators needs to be 2."))
    end
end


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
    h, σ = get_system_operators(N)

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
                    order=order)
    # find states which are on the equation rhs
    missing_states = find_missing(eqs)
    # set them to 0
    missing_subs = Dict(missing_states .=> 0)

    eqs = simplify(substitute(eqs, missing_subs))


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

end
