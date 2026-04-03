using Cairo, Compose
include("Lattice.jl")
include("QMC.jl")
include("QMC_plot.jl")

##v3 works but the lines are a little off

# triangle graph with defects
h = 4 #8  #use even numbers to avoid edges sticking out.
l = 4 #12
n = h*l

g = TriangleGrid(l,h)
A = Graphs.adjacency_matrix(g);


# approximate max cut

val2, X2 = QMC_SOC(A,     pauli_NPA1=true, verbose=false);

#print(X2)

X2[12, 16]=-0.75
X2[7, 8]=-0.75
X2[8, 12]=-0.75
X2[12, 15]=-0.75

for i=1:16
	for j=1:16
		X2[i, j]=-0.75

	end
end

for i=1:16 
	X2[i, i]=1
end

obj, l_edges=rounded_obj(X2, A)

print(l_edges)

plot2 = QMC_Plot(g,X2, large_edges=l_edges)
draw(PDF(join([n,"my_fig.pdf"]), 16cm, 16cm), plot2)
