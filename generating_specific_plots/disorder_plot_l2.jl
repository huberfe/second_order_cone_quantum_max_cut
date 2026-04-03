include("../QMC.jl")
include("../QMC_plot.jl")
include("../Lattice.jl")
include("../SwapMat.jl")
using Plots
using LaTeXStrings
using LinearAlgebra
using SparseArrays
using Arpack


function make_disorder_plot()
	l = h = 2

	J_diag = 1
	# open boundary lattice only for plotting
	A_open = ShastrySutherland2D(l,h; J=1, J_diag=1, periodic=false); 
	size(A_open)

	vals=zeros(length(LinRange(-3, -0.3, 20)))
	SOC_vals=zeros(length(LinRange(-3, -0.3, 20)))
	SOC_P1_vals=zeros(length(LinRange(-3, -0.3, 20)))
	c=0

	for log_sig in LinRange(-3, -0.3, 20)
	    println(log_sig)
	    A = ShastrySutherland2D_disorder(l,h; J=J, J_diag=J_diag, periodic=true, sig=10^log_sig);
	    val, X = QMC_SOC(A, verbose=false, pauli_NPA1=false, solver=solver)
	    c+=1
	    vals[c]=J
	    SOC_vals[c]=val
	    val, X = QMC_SOC(A, verbose=false, pauli_NPA1=true, solver=solver, sig=10^log_sig);)
	    SOC_P1_vals[c]=val
	    exact_vals[c]=-Arpack.eigs(-HamiltonianMaxCut(A, is_sp=true), nev=1, which=:LM)[1][1]
	end
	pl = Plots.plot(xlabel=L"$J/J_D$", ylabel=L"$\langle H \rangle$", title="Optimal Value of Relaxations");

    Plots.plot!(pl,vals,SOC_vals, label="SOC");
    Plots.plot!(pl,vals,SOC_P1_vals, label="SOC+P1");
    Plots.plot!(pl,vals,exact_vals, label="Exact Diag.");

    #Plots.display(pl)

    Plots.savefig(pl,"/ascldap/users/kevthom/quantum-max-cut/code/figs/disorder_plot_l2.png")


end
