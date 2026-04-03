include("QMC.jl")
include("QMC_plot.jl")
include("Lattice.jl")
using Plots
using LaTeXStrings
using LinearAlgebra
using FileIO, JLD2, Arpack


function opt_SS()
    ##This is for plotting the energy you get as a function of J

    l = h = 2

    J_diag = 1
    # open boundary lattice only for plotting
    A_open = ShastrySutherland2D(l,h; J=1, J_diag=1, periodic=false); 


    vals=zeros(length(collect(0.3:0.05:0.8)))
    SOC_vals=zeros(length(collect(0.3:0.05:0.8)))
    SOC_P1_vals=zeros(length(collect(0.3:0.05:0.8)))
    exact_vals=zeros(length(collect(0.3:0.05:0.8)))
    rounded_vals=zeros(length(collect(0.3:0.05:0.8))) 
    swap_vals=zeros(length(collect(0.3:0.05:0.8))) 
    c=0

    for J in collect(0.3:0.05:0.8)
        println(J)
        A = ShastrySutherland2D(l,h; J=J, J_diag=J_diag, periodic=true);
        #val, X = QMC_SOC(A, verbose=false, pauli_NPA1=false, solver="mosek")
        #val, X = QMC_SOC(A, verbose=false, pauli_NPA1=true, solver="mosek")

        c+=1
        vals[c]=J
        SOC_vals[c]=QMC_SOC(A, verbose=false, pauli_NPA1=false, solver="mosek")[1]
        val, X=QMC_SOC(A, verbose=false, pauli_NPA1=true, solver="mosek")
        SOC_P1_vals[c]=val 
        rounded_vals[c]=rounded_obj(X, A, samps=1000, var_bench_norm=false)
        exact_vals[c]=-Arpack.eigs(-HamiltonianMaxCut(A, is_sp=true), nev=1, which=:LM)[1][1]
        swap_vals[c]=l1_perm_program(A; solver="mosek", verbose=false, VarBench_norm=false, time_lim=2000)[1]
        #val=-Arpack.eigs(-HamiltonianMaxCut(A, is_sp=true), nev=1, which=:LM)[1][1]
        #val4, X = l1_perm_program(A; solver="mosek", verbose=false, VarBench_norm=false, time_lim=600)
        #swap_vals[c]=val4

    end
    dict1=Dict("vals"=> vals, "SOC_vals"=> SOC_vals, "SOC_P1_vals"=> SOC_P1_vals, "exact_vals"=> exact_vals, "rounded_vals"=> rounded_vals, "swap_vals"=> swap_vals)
    FileIO.save("SS_disc_vars_v2.jld2","dict1",dict1)
    FileIO.save("SS_disc_vars_v2_stash.jld2","dict1",dict1)

    pl = Plots.plot(xlabel=L"$J/J_D$", ylabel=L"$\langle H \rangle$", title="Optimal Value of Relaxations");

    Plots.plot!(pl,vals,SOC_vals, label="SOC");
    Plots.plot!(pl,vals,SOC_P1_vals, label="SOC+P1");
    Plots.plot!(pl, vals, rounded_vals, label="Rounded Objective Value")
    Plots.plot!(pl, vals, exact_vals, label="Exact Value")

    #Plots.display(pl)

    #Plots.savefig(pl,"/mnt/c/users/kevthom/Documents/gitlab_projects/quantum-max-cut/code/notebooks lattices/opt_SS_model.png")
    Plots.savefig(pl,"/ascldap/users/kevthom/quantum-max-cut/code/figs/ss_disc_w_rounded.pdf") 

    ##add optimal objective?  Will need to run on the big computer.
    #println(HamiltonianMaxCut(A))
end

#opt_SS()