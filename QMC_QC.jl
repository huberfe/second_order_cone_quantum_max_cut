using JuMP, LinearAlgebra, SparseArrays
using MosekTools
try import Mosek; catch err; println("MOSEK not installed"); end
try import SCS;   catch err; println("SCS not installed"); end
try import SDPA;  catch err; println("SDPA not installed"); end
try import Clarabel; catch err; println("Clarabel not installed"); end
import IterTools
include("SwapMat.jl")

### classical to quantum Max Cut
# NOT FINISHED

function QMCJ(A; J=1, pauli_NPA1=false, verbose=false, VarBench_norm=true)
    """ A: n x n adjacency matrix of Heisenberg Hamiltonian containing weights
    returns lower bound on energy of Hamiltonian H = sum_{ij} J*(XX + YY) + ZZ

    Natively in VARBENCH normalization!

    flags (true/false):
    pauli_NPA1 :  use Pauli NPA hierarchy level 1   # TODO!
    verbose    :  print solver info
    VarBench_norm  :  use normalization of Variational Benchmark paper
    H = sum_{ij in E} (XX + YY + ZZ)_ij
    https://arxiv.org/abs/2302.04919 / https://github.com/varbench/varbench/

    This optimization uses variable convention: X_ij = <(ij)>_rho
    not (yet) implemented: SWAP quantum Lasserre

    todo: dependency on weights, matrix A_xy and A_z
    """

    n  = size(A)[1]
    nE = sum(A) / 2 # number edges
    model = Model(JuMP.optimizer_with_attributes(Mosek.Optimizer))

    if verbose==false
        set_silent(model)
    end

    @variable(model, -1 <= Q[i=1:n, j=i+1:n] <= 1)      # variables X_ij with i<j
    @variable(model, -1 <= C[i=1:n, j=i+1:n] <= 1)      # variables C_ij with i<j


    ##### SDP ####
    # sym red 3-RDM constraint
    for (i,j,k) in IterTools.subsets(1:n, 3) #ordered subsets i<j<k
        a = 1 + C[i,j] + C[j,k] + C[i,k]
        b = 1 + C[i,j] - C[j,k] - C[i,k]
        c = 1 - C[i,j] + C[j,k] - C[i,k]
        d = 1 - C[i,j] - C[j,k] - C[i,k]
        ee= Q[i,j]
        f = Q[i,k]
        g = Q[j,k]

        @constraint(model, 0 ≤ a)
        @constraint(model, a + b + c + d == 1)

        W = [b  ee f;
             ee c  g;
             f  g  d]
        @constraint(model,  W in PSDCone())
    end

    #TODO: correct constraints!
    ##### Pauli NPA ####
    if pauli_NPA1 == true
        @variable(model, M[1:n,1:n] in PSDCone())      # variables X_ij with i<j Symmetric
        for i=1:n
            @constraint(model, M[i,i] == 3)  #diag = <X^2 + Y^2 + Z^2> = <3I> = 3
        end
        for i=1:n, j=i+1:n
            @constraint(model, M[i,j] == 2*X[i,j] - 1 )  # <XX> + <YY> + <ZZ> = 2*swap - II
        end
    end

    # translate to swap normalization:
    for i=1:n, j=i+1:n
        X[i,j] = (1 + C[i,j] + Q[i,j]) / 2   # swap = 1/2(II + XX + YY + ZZ)
    end

    #### Energy ####
    #only access i<j adjacency matrix part
    # - sum_ij tr[(1 - swap_ij)/2 rho] = - sum_ij tr[singlet_ij rho]

    if VarBench_norm == false
        H = - (nE - sum([X[i,j]*A[i,j] for i=1:n, j=1:n if i<j])) / 2
            elseif VarBench_norm == true
            #H = 2 * sum([X[i,j]*A[i,j] for i=1:n, j=1:n if i<j]) - nE
        end

        @objective(model, Min, H)
        optimize!(model)

        if verbose == true
            #print(model);
            @show(solution_summary(model)); @show(termination_status(model)); @show(objective_value(model))
        end
        return objective_value(model), JuMP.value.(C), JuMP.value.(Q), JuMP.value.(X)
        end
