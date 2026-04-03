##really big SS with disorder

using Cairo, Compose, FileIO, JLD2, LinearAlgebra, LaTeXStrings
include("../Lattice.jl")
include("../QMC.jl")
include("../QMC_plot.jl")
include("../SwapMat.jl")

function make_SS_plot()
	##this is the case that will make a 256 qubit instance.
	l=h=8
	J_diag=1
	J=0.4
	A=ShastrySutherland2D_disorder(l,h; J=J, J_diag=J_diag, periodic=true, sig=0.05)
	val, X = QMC_SOC(A, verbose=true, pauli_NPA1=true, solver="scs")
	obj, l_edges=rounded_obj(X, A)
	dict1=Dict("X"=> X, "val"=> val, "A"=> A, "l"=> l, "l_edges"=>l_edges)
	FileIO.save("big_SS_varsv2l8"*string(round(Int, time()))*".jld2","dict1",dict1)
	FileIO.save("big_SS_vars_stashv2l8"*string(round(Int, time()))*".jld2","dict1",dict1)

end

function mk_heat_map()
	samps=200

	dict2=FileIO.load("big_SS_varsv2l81771609772.jld2", "dict1")
	X=get(dict2, "X", "Sorry, no X")
	val=get(dict2, "val", "Sorry, val")
	A=get(dict2, "A", "Sorry, no A")
	l=get(dict2, "l", "Sorry, no l")
	l_edges=dict2["l_edges"]
	h=l
	A_open = shastry_sutherland2D(l,h, J=1, J_diag=1, periodic=false); 
	for i in 1:samps
		plot = QMC_Plot(A_open,X, large_edges=l_edges)
		draw(SVG("big_SS_with_disorder_heat_mapv3"*string(i)*".svg", 16cm, 16cm), plot)

		
	end


end


mk_heat_map()
