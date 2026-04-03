include("../QMC.jl")
include("../QMC_plot.jl")
include("../Lattice.jl")
include("../SwapMat.jl")
using FileIO, JLD2, Arpack, SparseArrays, LinearAlgebra, LaTeXStrings, Plots


function make_disorder_plot(;samps=100)
	l = h = 2
	solver="mosek"
	J=0.4
	num_points=5

	J_diag = 1
	# open boundary lattice only for plotting
	A_open = ShastrySutherland2D(l,h; J=1, J_diag=1, periodic=false); 
	size(A_open)

	vals=zeros(length(LinRange(-3, -0.3, num_points)))
	SOC_vals=zeros(length(LinRange(-3, -0.3, num_points)))
	SOC_P1_vals=zeros(length(LinRange(-3, -0.3, num_points)))
	swap_vals=zeros(length(LinRange(-3, -0.3, num_points)))
	c=0

	for log_sig in LinRange(-3, -0.3, num_points)
	    println(log_sig)

	    avg_SOC_dev=0
	    avg_SOC_P1_dev=0
	    avg_swap_dev=0
	    for b in collect(1:1:samps)
	    	if mod(b, 10)==0
				println(log_sig)
				println(b)
			end
	    	A = ShastrySutherland2D_disorder(l,h; J=J, J_diag=J_diag, periodic=true, sig=10^log_sig);
	    	val1, X = QMC_SOC(A, verbose=false, pauli_NPA1=false, solver=solver)
	    	val2, X = QMC_SOC(A, verbose=false, pauli_NPA1=true, solver=solver)
	    	val4, X = l1_perm_program(A, time_lim=10000)
	    	val3= -Arpack.eigs(-HamiltonianMaxCut(A, is_sp=true), nev=1, which=:LM)[1][1]
	    	avg_SOC_dev+= abs(val1-val3)/(abs(val3)*samps)
			avg_SOC_P1_dev+= abs(val2-val3)/(abs(val3)*samps)
			avg_swap_dev+= abs(val4-val3)/(abs(val4)*samps)


	    end
	    c+=1
	    vals[c]=10^log_sig
	    SOC_vals[c]=avg_SOC_dev;
		SOC_P1_vals[c]=avg_SOC_P1_dev;
		swap_vals[c]=avg_swap_dev

	end
	dict1=Dict("vals"=> vals, "SOC_vals"=> SOC_vals, "SOC_P1_vals"=> SOC_P1_vals, "swap_vals"=> swap_vals)
	FileIO.save("disorder_vars.jld2","dict1",dict1)
	FileIO.save("disorder_vars_stash.jld2","dict1",dict1)

	pl = Plots.plot(xlabel=L"$\sigma$",  ylabel=L"\mathbb{E}|\langle H\rangle -\lambda_{min}(H)|/|\lambda_{min}(H)|", title="Expected Discrepancy with Disorder",  xaxis=:log);

    Plots.plot!(pl,vals,SOC_vals, label="SOC");
    Plots.plot!(pl,vals,SOC_P1_vals, label="SOC+P1");
    Plots.plot!(pl,vals,swap_vals, label="Lvl-1 SWAP");

    #Plots.display(pl)

    Plots.savefig(pl,"/ascldap/users/kevthom/quantum-max-cut/code/figs/disorder_plot_l2_disc.pdf")
    #Plots.savefig(pl,"/mnt/c/users/kevthom/Documents/gitlab_projects/quantum-max-cut/code/disorder_plot_l2_disc.pdf")


end

make_disorder_plot(samps=1000)
