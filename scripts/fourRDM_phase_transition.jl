include("../QMC.jl")
include("../QMC_plot.jl")
include("../Lattice.jl")
include("../SwapMat.jl")
using FileIO, JLD2, Plots, LinearAlgebra, LaTeXStrings
import Cairo # Cairo is needed for Compose.PDF
import Compose #import Cairo


##goal is to investigate if 4RDMs change the phase transition.

function make_bunch_of_heat_maps(J_vals)
	l = h = 4
	# open boundary lattice only for plotting
		A_open = ShastrySutherland2D(l,h; J=1, J_diag=1, periodic=false); 
		solver="mosek"
		FS = ShastrySutherland2D_4RDM(l,h)


	for J in J_vals 
		println(J)
		A = ShastrySutherland2D(l,h; J=J, J_diag=1, periodic=true);
		val, X = QMC_SOC(A, verbose=true, pauli_NPA1=false, four_qb = FS, solver=solver)
		dict1=Dict("val"=> val, "X"=> X)
		FileIO.save("SS_marg_vars" * string(J) * ".jld2","dict1",dict1)
		FileIO.save("copy_SS_marg_vars" * string(J) * ".jld2","dict1",dict1)
		plot = QMC_Plot(A_open,X)
		name = "ShastrySutherland_marg_" * string(l) * string(h) * "_" * string(J) * ".pdf"
		Compose.draw(Compose.PDF(name, 16Compose.cm, 16Compose.cm), plot) #print to file
	end

end

make_bunch_of_heat_maps([0.45, 0.47, 0.5, 0.52, 0.55])
