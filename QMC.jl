using JuMP, LinearAlgebra, SparseArrays
using MosekTools
try import Mosek; catch err; println("MOSEK not installed"); end
try import SCS;   catch err; println("SCS not installed"); end
try import SDPA;  catch err; println("SDPA not installed"); end
try import Clarabel; catch err; println("Clarabel not installed"); end
import IterTools
include("SwapMat.jl")

function QMC_SOC(A; pauli_NPA1=false, four_qb=false, partition_31=true, partition_31_relax=false,
                    starbound=false, starbound4=false, SDP_form=false,
                    solver="mosek", verbose=false, VarBench_norm=false)
    
    """ Given an adjacency matrix A, QMC_SOC provides a lower bound on the ground state energy
        of the Quantum Max Cut Hamiltonian:
        H = - sum_{ij in E} 1/2 * (1 - swap_ij) ,     where swap = 1/2*(II + XX + YY + ZZ) 
          = - sum_{ij in E} 1/4 * (1 - XX - YY - ZZ)

        corresponds to
        H = - sum_{ij in E} 1/2 * (1 - x_ij)
        > x_ij = - 1: singlet
        > x_ij = + 1: orthogonal to singlet  (Opposite to Robbie King notation)

        It returns the lower bound to min tr(H \rho) and 
        a set of (non-physical) edge weights approximating X_ij = <(ij)>_rho,
        obtained by optimizing over mutually compatible 3-body reduced density matrices.

        Arguments:
        A: n x n adjacency matrix of Heisenberg Hamiltonian containing weights
        pauli_NPA1 :  bool; use Pauli NPA hierarchy level 1
        four_qb :     use 4 qubit constraints. true, false, or a list of (ordered) lists 
                      with 4 elements each. E.g. [[0,2,3,5], [2,3,5,6]
        partition_31 : bool; use partition [3,1] constraint of the 4-qb constraint
        partition_31_relax : bool;  true: uses SOC relaxation on 2x2 minors
                                    false: use SDP form of [3,1] partition 3x3 matrix constraint. 
                             
        starbound  :  bool; use the star bound
        starbound4 :  bool; use Thompson starbound on 4 RDM
        SDP_form   :  bool; use 2x2 SDP formulation of Second order cone for 3RDM.
        verbose    :  print solver info
        VarBench_norm  :  use normalization of Variational Benchmark paper 
                          H = sum_{ij in E} (XX + YY + ZZ)_ij
                          https://arxiv.org/abs/2302.04919 / https://github.com/varbench/varbench/
        solvers : "mosek", "scs", "sdpa", "clarabel"
            mosek: https://www.mosek.com/
            SCS: https://www.cvxgrp.org/scs/
            SDPA: https://github.com/jump-dev/SDPA.jl
            ECOS: https://github.com/embotech/ecos      (only SOC!)
            Clarabel: https://clarabel.org/stable/julia/jump/

        This optimization uses variable convention: X_ij = <(ij)>_rho
        where (ij) is the swap operator
        Returns:
        lower bound : float
        edge values:  dict of floats
    """
    
    n  = size(A)[1]
    nE = sum(A) / 2 # number of (weighted) edges

    # select solver
    if solver=="mosek"
        model = Model(JuMP.optimizer_with_attributes(Mosek.Optimizer))
    elseif solver=="scs"
        model = Model(SCS.Optimizer)
    elseif solver=="sdpa"
        model = Model(SDPA.Optimizer)
        set_attribute(model, "Mode", SDPA.PARAMETER_DEFAULT)
    elseif solver=="clarabel"
        model = JuMP.Model(Clarabel.Optimizer)
    end
    
    if verbose==false
        set_silent(model)
    end
    
    @variable(model, -1 <= X[i=1:n, j=i+1:n] <= 1)      # variables X_ij with i<j

    ##### Second order cone #####
    # sym red 3-RDM constraint
    for (i,j,k) in IterTools.subsets(1:n, 3) #ordered subsets i<j<k
        rp =  (X[i,j] + X[j,k] + X[i,k]) / 3
        r0 =  (3 - X[i,j] - X[j,k] - X[i,k]) / 3
        r1 =  (2*X[j,k] - X[i,k] - X[i,j]) / 3
        r2 =  (X[i,j] - X[i,k]) / sqrt(3)

        @constraint(model, 0 ≤ rp ≤ 1)
        @constraint(model, 0 ≤ r0 ≤ 1)
        @constraint(model, rp + r0 == 1)

        if SDP_form == false
            @constraint(model, [r0,r1,r2] in SecondOrderCone())
        else
            W = [(r0 + r1)  r2; 
                 r2  (r0 - r1)]
            @constraint(model,  W in PSDCone())    
        end
    end


    ##### Four qubit positivity constraints ####
    if four_qb != false
        if four_qb == true 
            L = IterTools.subsets(1:n,4)    # all involutions of form (ij)(kl) = all ordered subsets i<j<k<l
        else
            L = four_qb  # use provided list of 4-RDMs
        end
        
        @variable(model, -1 <= Y[i=1:n, j=1:n, k=1:n, l=1:n] <= 1)      
        # Y[i,j,k,l] stand for (ij)(kl) expectation values, in that order
        # condition: i<j, k<l. Note: too many vars defined
    
        # impose 4-RDM constraints on subset ijkl
        for (i,j,k,l) in L
            
            #partition [4]:
            R4 = 1/12*(-3 + 2*(X[i,j] + X[i,k] + X[i,l] + X[j,k] + X[j,l] + X[k,l]) 
                          + Y[i,j,k,l] + Y[i,k,j,l] + Y[i,l,j,k])
            
            #partition [3,1]:
            R31 = 1/4*( 3 - Y[i,j,k,l] - Y[i,k,j,l] - Y[i,l,j,k])
            
            #partition [2,2]:
            R22 =  1/6*( 3 - X[i,j] - X[i,k] - X[i,l] - X[j,k]  - X[j,l]  - X[k,l]
                    + Y[i,j,k,l]  + Y[i,k,j,l] + Y[i,l,j,k] )
                
            B00 = (3 + X[i,j] - 2 * X[i,k] - 2 * X[i,l] - 2 * X[j,k] - 2 * X[j,l] + X[k,l] 
                 - Y[i,j,k,l] + 2 * Y[i,k,j,l] + 2 * Y[i,l,j,k])
            B01 = sqrt(3)*( - X[i,k] + X[i,l] + X[j,k] - X[j,l] 
                 + Y[i,k,j,l] - Y[i,l,j,k] )
            B11 = 3*(1 - X[i,j] - X[k,l] 
                 + Y[i,j,k,l] )
            
            BB  = [B00 B01;
                   B01 B11]  # matrix SDP form

            b0 = (B00 + B11)/2
            b1 = (B00 - B11)/2
            b2 = B01         # to formulate SOC cone

            # positive and normalized irreps projectors
            @constraint(model, 0 ≤ R4  ≤ 1)  #irrep [4]
            @constraint(model, 0 ≤ R31 ≤ 1)
            @constraint(model, 0 ≤ R22 ≤ 1)
            @constraint(model, R4 + R31 + R22 == 1)

            #positivity on [2,2] and [3,1] irreps
            # irrep [2,2]
            if SDP_form == true
                @constraint(model, BB in PSDCone())   
            else
                @constraint(model, [b0,b1,b2] in SecondOrderCone())   # 2x2 matrix B is a SOC
            end

            if partition_31 == true
                A00 = 2/3*( 3 + 2*( X[i,j] + X[i,k] - X[i,l] + X[j,k] - X[j,l] - X[k,l])  
                    - Y[i,j,k,l] - Y[i,k,j,l] - Y[i,l,j,k] )
                A01 = sqrt(2)/3*( -2*X[i,j] + X[i,k] - X[i,l] + X[j,k] - X[j,l] + 2 *X[k,l] 
                     + 4 * Y[i,j,k,l]  - 2 * Y[i,k,j,l] - 2 * Y[i,l,j,k] )
                A02 = sqrt(6)/3*( - X[i,k] - X[i,l] + X[j,k] + X[j,l] 
                     + 2 * Y[i,k,j,l] - 2 * Y[i,l,j,k] )
                A11 = 2/3*( 3 + X[i,j] - X[k,l] + 2 * X[j,l] - 2 * X[i,k] + 2 * X[i,l] - 2 * X[j,k] 
                     + Y[i,j,k,l] - 2 * Y[i,k,j,l] - 2 * Y[i,l,j,k] )
                A12 = 2/sqrt(3)*( -X[i,k] - X[i,l] + X[j,k] + X[j,l] 
                     - Y[i,k,j,l] + Y[i,l,j,k] )
                A22 = 2*( 1 - X[i,j] - Y[i,j,k,l] + X[k,l] )
                
                AA  = [A00 A01 A02; 
                       A01 A11 A12; 
                       A02 A12 A22]
            
                if partition_31_relax == false
                    # no relaxation, use SDP constraint
                    @constraint(model, AA in PSDCone())
                else
                    # SOC relaxation of AA positive semidefinite:
                    #AA  = [A00 A01 A02; 
                    #       A01 A11 A12; 
                    #       A02 A12 A22]
                    # 01-block
                    AA_a0 = (A00 + A11)/2
                    AA_a1 = (A00 - A11)/2
                    AA_a2 = A01
    
                    # 02-block
                    AA_b0 = (A00 + A22)/2
                    AA_b1 = (A00 - A22)/2
                    AA_b2 = A02
    
                    # 12-block
                    AA_c0 = (A11 + A22)/2
                    AA_c1 = (A11 - A22)/2
                    AA_c2 = A12
    
                    # SOC relaxation of psd-ness of AA: all 2x2 minors psd
                    @constraint(model, 0 ≤ A00)
                    @constraint(model, 0 ≤ A11)
                    #@constraint(model, 0 ≤ A22)
                    @constraint(model, [AA_a0, AA_a1, AA_a2] in SecondOrderCone())
                    @constraint(model, [AA_b0, AA_b1, AA_b2] in SecondOrderCone())
                    @constraint(model, [AA_c0, AA_c1, AA_c2] in SecondOrderCone())
                end
            end
        end
    end

    ##### Pauli NPA ####
    if pauli_NPA1 == true
        @variable(model, M[1:n,1:n] in PSDCone())        # variables X_ij with i<j Symmetric
        for i=1:n
            @constraint(model, M[i,i] == 3)              #diag = <X^2 + Y^2 + Z^2> = <3I> = 3
        end
        for i=1:n, j=i+1:n
            @constraint(model, M[i,j] == 2*X[i,j] - 1 )  # <XX> + <YY> + <ZZ> = 2*swap - II
        end
    end


    # star bound : - sum_{j \in N(i)} X_{ij} \leq 1
    # careful: opposite notation than in R. King, , https://arxiv.org/abs/2209.02589, page 7 Lemma 7.
    for i = 1:n
        if starbound == true
            @constraint(model, - sum([X[i,j]*A[i,j] for j=1:n]) ≤ 1)
        end
    end

    # star bound by Thompson on 4 RDM. TODO: ADD reference
    if starbound4 == true
        for (i,j,k,l) in IterTools.subsets(1:n, 4) #ordered subsets i<j<k<l
            # X_{12} + X_{13}+X_{14} >= -2/3
            @constraint(model, (X[i,j] + X[i,k] + X[i,l]) ≥ -2/3)
        end
    end
    

    #### Energy ####
    #only access i<j adjacency matrix part
    # H = - sum_ij tr[(1 - swap_ij)/2 rho]
    # X_ij = <(ij)>_rho 
    #      = 1/2*(1 + <XX> + <YY> + <ZZ>)

    if VarBench_norm == false
        H = - (nE - sum([X[i,j]*A[i,j] for i=1:n, j=1:n if i<j])) / 2
    elseif VarBench_norm == true
        H = 2 * sum([X[i,j]*A[i,j] for i=1:n, j=1:n if i<j]) - nE
    end

    @objective(model, Min, H) 
    optimize!(model)
    
    if verbose == true
        #print(model);  
        @show(solution_summary(model)); @show(termination_status(model)); @show(objective_value(model))
    end
    return objective_value(model), JuMP.value.(X)
end



function HamiltonianMaxCut(A; is_sp=false::Bool)
    """ exact Ground state by diagonalization """
    if is_sp
        n = size(A)[1]
        HH = SparseArrays.spzeros(2^n, 2^n) 
        iden=SparseArrays.sparse(LinearAlgebra.I, 2^n, 2^n)
        for i=1:n, j=i+1:n
            if abs(A[i,j])>10^(-5)
                HH = HH - A[i,j]*(iden - swap_padded(i,j,n, is_sp=is_sp))/2
            end
        end
        return HH
    else  
        # exact Ground state by diagonalization
        n = size(A)[1]
        HH = 0*LinearAlgebra.I
        for i=1:n, j=i+1:n
            if abs(A[i,j])>10^(-5)
                HH = HH - A[i,j]*(LinearAlgebra.I - swap_padded(i,j,n))/2
            end
        end
        return HH


    end
    # exact Ground state by diagonalization
    #n = size(A)[1]
    #HH = 0*I
    #for i=1:n, j=i+1:n
    #    if A[i,j] == 1
    #        HH = HH - (I - swap_padded(i,j,n))/2
    #    end
    #end
    #return HH
end


function QMC_to_varbench(A, val)
    """ transforms QMC to Variational Benchmark normalization """
    nE = sum(A) / 2
    return (val + nE/4)*4
end





function rounded_obj(X, A; samps=1000, var_bench_norm=true)
    """ X a matrix of size binom(n,2) X[i,j] is the edge value returned from  SOC + Lasserre-1. 
        i.e. is the output of QMC_SOC
        A is the matrix corresponding to the problem description, 
        returns avg_obj/samps, large_edges 
    """

    ##first convert X to dense
    n=size(A)[1]
    t=0.771
    nE = sum(A) / 2
    X_dense=zeros(n, n)
    for key in eachindex(X)
        X_dense[key[1], key[2]]=X[key]
    end
    X=X_dense+X_dense'
    ##now need to construct moment matrix
    M=zeros(n, n)
    for i=1:n, j=i+1:n
        M[i, j]=2*X[i, j]-1
    end
    M=M+M'+3*LinearAlgebra.I
    ##need Cholesky decomposition.  First compute eigendecomposition of M then throw away negative eigenvalues and use decomposition
    ##to make Cholesky vectors
    vals, vecs=LinearAlgebra.eigen(M)
    root_vals=zeros(n)
    for j=1:length(vals)
        if vals[j] < 0
            root_vals[j]=0
        else 
            root_vals[j]=sqrt(vals[j])
        end
    end
    V= vecs*diagm(root_vals)
    
    ##find edges with large SDP values. I dont think I need the actual edges just the total weight
    large_vertices=[]
    large_edges=[]
    big_W=0
    for i=1:n, j=i+1:n 
        if abs(A[i, j])>= 10^(-5)
            if 1-2t > X[i, j]
                push!(large_edges, (i, j))
                big_W+=A[i, j]

                if !(i in large_vertices)
                    push!(large_vertices, i)
                end 
                if !(j in large_vertices)
                    push!(large_vertices, j)
                end 

            end
        end
    end
    
    ##edges adjacent to large edges
    ad_W=0
    for i=1:n, j=i+1:n 
        if abs(A[i, j])>= 10^(-5)
            if (i in large_vertices) && !(j in large_vertices)
                ad_W+=A[i, j]
            elseif !(i in large_vertices) && (j in large_vertices)
                ad_W+=A[i, j]
            end 
        end 
    end 
    A_small=zeros(n, n)
    small_W=0
    for i=1:n, j=i+1:n 
        if !(i in large_vertices) && !(j in large_vertices)
            A_small[i, j]=A[i, j] 
            small_W+=A[i, j]
        end 
    end 
    A_small=A_small+A_small'

    avg_obj=0
    for j=1:samps 
        ##vecs*diagm(root_vals)*(vecs*diagm(root_vals))^T =M.  The vectors we want are the rows of vec*diagm(root_vals)
        R=randn(n, 3)
        
        #println(V)
        #println(V*V')
        #println(M)
        inv_norms=1 ./ sqrt.(LinearAlgebra.diag(V*R*R'*V'))
        bloch=diagm(inv_norms)*V*R
        bloch_sq=bloch * bloch'
        if var_bench_norm
            cur_prod_obj=LinearAlgebra.tr(A * bloch_sq)/2
            cur_match_obj= -3*big_W+LinearAlgebra.tr(A_small' * bloch_sq)/2
            
        else 
            ##########################still need to code this part
            cur_prod_obj= -nE/4+LinearAlgebra.tr(A * bloch_sq)/8
            cur_match_obj= -big_W-ad_W/4-small_W/4+LinearAlgebra.tr(A_small' * bloch_sq)/8

        end 
        
        if cur_prod_obj <= cur_match_obj
            avg_obj+=cur_prod_obj
        else 
            avg_obj+=cur_match_obj
        end 
    end
    #println(avg_obj/samps)
    return avg_obj/samps, large_edges
end

#A=[0 1 1 1 1;1 0 1 1 1; 1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0]
#A=[0 1 0; 1 0 0;0 0 0]
#X=QMC_SOC(A, solver="scs", pauli_NPA1=true)[2]
#print(rounded_obj(X, A))

function l1_perm_program(A; solver="mosek", verbose=false, VarBench_norm=false, time_lim=600)
    n  = size(A)[1]
    nE = sum(A) / 2 # number edges
    n_pairs=Int(n*(n-1)/2)
    edge_list=collect(IterTools.subsets(1:n, 2))
    #println(edge_list)
    ##first index will correspond to I, rest of them will be edges.

    # select solver
    if solver=="mosek"
        model = Model(JuMP.optimizer_with_attributes(Mosek.Optimizer))
    elseif solver=="scs"
        model = Model(SCS.Optimizer)
    #elseif solver=="ecos"
    #    model = Model(ECOS.Optimizer)
    #    set_attribute(model, "maxit", 100)
    elseif solver=="sdpa"
        model = Model(SDPA.Optimizer)
        set_attribute(model, "Mode", SDPA.PARAMETER_DEFAULT)
    elseif solver=="clarabel"
        model = JuMP.Model(Clarabel.Optimizer)
    end

    if verbose==false
        set_silent(model)
    end

    set_time_limit_sec(model, time_lim)

    @variable(model, X[1:n_pairs+1,1:n_pairs+1] in PSDCone())
    @constraint(model, X[1,1] == 1 )
    for i=2:n_pairs+1
        @constraint(model, X[i,i] == 1 )
        for j=i+1:n_pairs+1
            cond1=(edge_list[i-1][1] in edge_list[j-1]) && !(edge_list[i-1][2] in edge_list[j-1])
            cond2=!(edge_list[i-1][1] in edge_list[j-1]) && (edge_list[i-1][2] in edge_list[j-1])
            if cond1 || cond2
                
                inds=sort(union(edge_list[i-1], edge_list[j-1]))
                #println(inds)
                
                e1=[inds[1], inds[2]]
                e2=[inds[2], inds[3]]
                e3=[inds[1], inds[3]]
                ind1=findall(x->x==e1, edge_list)
                ind2=findall(x->x==e2, edge_list)
                ind3=findall(x->x==e3, edge_list)
                #println(e1)
                #println(ind1)
                #println(e2)
                #println(ind2)
                #println(e3)
                #println(ind3)
                @constraint(model, X[i,j] == (X[1, ind1[1]+1]+X[1, ind2[1]+1]+X[1, ind3[1]+1]-1)/2 )
        
            end 
        end
    end

    ##now for the objective function
    H=0
    for i=1:n, j=i+1:n 
        if abs(A[i, j])>= 10^(-5)
            ind=findall(x->x==[i, j], edge_list)
            if VarBench_norm == false
                H+=-A[i, j]*(1-X[1, ind[1]+1])/2
            else 
                H+= (2*X[1, ind[1]+1]-1)*A[i, j]
            end
        end
    end 

    @objective(model, Min, H) 
    optimize!(model)
    #println(solve_time(model))
    
    if verbose == true
        #print(model);  
        @show(solution_summary(model)); @show(termination_status(model)); @show(objective_value(model))
    end
    return objective_value(model), JuMP.value.(X)

end 

#A=[0 1 1 1 1;1 0 1 1 1; 1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0]
#A=[0 1 1 1;1 0 1 1; 1 1 0 1; 1 1 1 0]
#A=[0 1 0.1; 1 0 0; 0.1 0 0]
#println(l1_perm_program(A, VarBench_norm=true))
