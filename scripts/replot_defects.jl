using Cairo, Compose, FileIO, JLD2
include("../Lattice.jl")
include("../QMC.jl")
include("../QMC_plot.jl")

dict1=FileIO.load("defect_vars1771531439461584e9.jld2", "dict1")
l=12
h=12
n=l*h

g = TriangleGrid(l,h)
A = Graphs.adjacency_matrix(g);


X2=dict1["X2"]
X4=dict1["X4"]
l_edges2=dict1["l_edges2"]
l_edges4=dict1["l_edges4"]
A=dict1["A"]
g=dict1["g"]
A_def=dict1["A_def"]
g_def=dict1["g_def"]




samps=100

for j=1:1:samps

	
	plot2 = QMC_Plot(g,X2, large_edges=l_edges2)
	draw(SVG(join([n,"triangle_SOC_Pauli"*string(j)*".svg"]), 16cm, 16cm), plot2)
end
for j=1:1:samps
	
	plot4 = QMC_Plot(g_def,X4, large_edges=l_edges4)
	draw(SVG(join([n,"triangle_def_SOC_Pauli"*string(j)*".svg"]), 16cm, 16cm), plot4)

end
