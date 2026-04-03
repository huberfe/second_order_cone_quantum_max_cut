include("../QMC.jl")
include("../QMC_plot.jl")
include("../Lattice.jl")
include("../SwapMat.jl")
using Plots
using LaTeXStrings
using LinearAlgebra
using SparseArrays
using Arpack
using FileIO, JLD2

function calc_er_graph(;samps=1000)
	ps=[0.2, 0.3, 0.5]
	max_n=14
	pl = Plots.plot(xlabel="n", ylabel=L"\mathbb{E}|\langle H\rangle -\lambda_{min}(H)|/|\lambda_{min}(H)|", title="Average Discrepancy", legend=:topleft);
	vals=collect(4:1:max_n)
	data=Dict{String,Any}("vals"=> vals)

	for p in ps	

		
		SOC_vals= zeros(length(collect(4:1:max_n)))
		SOC_P1_vals= zeros(length(collect(4:1:max_n)))
		swap_vals= zeros(length(collect(4:1:max_n)))
		c=0
		
		for n in collect(4:1:max_n)
			println(n)
			avg_SOC_dev=0
			avg_SOC_P1_dev=0
			avg_swap_dev=0
			for q in collect(1:1:samps)
				if mod(q, 100)==0

					println(q)
				end
				A=Erdos_renyi(n, p)
				val1, X = QMC_SOC(A, verbose=false, pauli_NPA1=false, solver="mosek")
				val2, X = QMC_SOC(A, verbose=false, pauli_NPA1=true, solver="mosek")
				val4, X = l1_perm_program(A; solver="mosek", verbose=false, VarBench_norm=false, time_lim=600)
				val3 = -Arpack.eigs(-HamiltonianMaxCut(A, is_sp=true), nev=1, which=:LM)[1][1]
				avg_SOC_dev+= abs(val1-val3)/(abs(val3)*samps)
				avg_SOC_P1_dev+= abs(val2-val3)/(abs(val3)*samps)
				avg_swap_dev+= abs(val4-val3)/(abs(val3)*samps)
			end
			c+=1
			vals[c]=n;
			SOC_vals[c]=avg_SOC_dev;
			SOC_P1_vals[c]=avg_SOC_P1_dev;
			swap_vals[c]=avg_swap_dev;



		end
		temp_dict=Dict{String,Any}("SOC_vals" * string(p)=> SOC_vals, "SOC_P1_vals" * string(p)=> SOC_P1_vals, "swap_vals" * string(p)=> swap_vals)
		merge!(data, temp_dict)


		Plots.plot!(pl,vals,SOC_vals, label="SOC, p="*string(p));
		Plots.plot!(pl,vals,SOC_P1_vals, label="SOC+P1, p="*string(p));
		Plots.plot!(pl,vals,swap_vals, label="Lvl-1 SWAP, p="*string(p));
		
	end
	FileIO.save("data_er_graph.jld2","data",data)
	FileIO.save("copy_data_er_graph.jld2","data",data)
	Plots.savefig(pl,"/ascldap/users/kevthom/quantum-max-cut/code/figs/opt_ER_point2_disc_v3.pdf")	
	#Plots.savefig(pl,"/mnt/c/users/kevthom/Documents/gitlab_projects/quantum-max-cut/code/opt_ER_point2_disc_v2.pdf")

end 

calc_er_graph(samps=5000)

