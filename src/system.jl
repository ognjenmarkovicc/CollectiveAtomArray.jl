"""
Code for setting up the system operators,
    and setting initial conditions.
"""

"""
    get_system_hilbert(N::Integer)
    
Get the Hilbert space for a system of N two-level atoms.
"""
function get_system_hilbert(N::Integer)
    # Hilbert space of the N two-level atoms
    # in this case NLevelSpace is a two-level space
    h = ⊗([NLevelSpace(Symbol(:atom, i), (:g, :e)) for i=1:N]...)

    return h
end

"""
    get_system_sigma(N::Integer)
    
Get the transition (sigma) operator
for a system of N two-level atoms.
"""
function get_system_sigma(N::Integer)

    h = get_system_hilbert(N)
    # operator for the transition from level j to level i for atom k
    # defined as a function
    σ(i,j,k) = Transition(h,Symbol(:σ_,k),i,j,k)

    return σ
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

    σ = get_system_sigma(N)

    if order == 2
        relevant_ops = [σ(:e, :g, i)*σ(:g, :e, j) for i in 1:N for j in 1:N if j>=i]
        append!(relevant_ops, [σ(:e, :e, i)*σ(:e, :e, j) for i in 1:N for j in 1:N if j>i])

        return relevant_ops 
    else
        throw(ArgumentError("order in get_nonzero_operators needs to be 2."))
    end
end