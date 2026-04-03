using Cairo, Compose, FileIO, JLD2
include("../Lattice.jl")
include("../QMC.jl")
include("../QMC_plot.jl")

# triangle graph with defects
h = 12 #8  #use even numbers to avoid edges sticking out.
l = 12 #12
n = h*l

g = TriangleGrid(l,h)
A = Graphs.adjacency_matrix(g);

#defect
g_def = DeleteRandomEdges(g,0.97)  #0.97
A_def = Graphs.adjacency_matrix(g_def);

#GraphPlot.gplot(g, layout=GraphPlot.spectral_layout)
#sprint layout (can be messed up too..)
#layout=(args...)->GraphPlot.spring_layout(args...; C=1)
#GraphPlot.gplot(g, layout=layout)

# approximate max cut
val1, X1 = QMC_SOC(A,     pauli_NPA1=false, verbose=true);
val2, X2 = QMC_SOC(A,     pauli_NPA1=true, verbose=true);

# with defects
val3, X3 = QMC_SOC(A_def, pauli_NPA1=false, verbose=true);
val4, X4 = QMC_SOC(A_def, pauli_NPA1=true, verbose=true);

dict1=Dict("X1"=> X1, "X2"=> X2, "X3"=> X3, "X4"=> X4)
FileIO.save("defect_vars.jld2","dict1",dict1)
FileIO.save("copy_defect_vars.jld2","dict1",dict1)

plot1 = QMC_Plot(g,X1)
draw(PDF(join([n,"triangle_SOC.pdf"]), 16cm, 16cm), plot1)

plot2 = QMC_Plot(g,X2)
draw(PDF(join([n,"triangle_SOC_Pauli.pdf"]), 16cm, 16cm), plot2)

plot3 = QMC_Plot(g_def,X3)
draw(PDF(join([n,"triangle_def_SOC.pdf"]), 16cm, 16cm), plot3)

plot4 = QMC_Plot(g_def,X4)
draw(PDF(join([n,"triangle_def_SOC_Pauli.pdf"]), 16cm, 16cm), plot4)
