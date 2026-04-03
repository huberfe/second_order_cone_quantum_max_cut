include("../QMC.jl")
include("../QMC_plot.jl")
include("../Lattice.jl")

import Cairo # Cairo is needed for Compose.PDF
import Compose #import Cairo
############################################

import Graphs
import GraphPlot
A = kagome(4,4, periodic=false);
g = Graphs.SimpleGraph(A)
layout=(args...)->GraphPlot.spring_layout(args...; C=1)
plot=GraphPlot.gplot(g, layout=layout)
Compose.draw(Compose.PDF("kagome_plot1.pdf", 16Compose.cm, 16Compose.cm), plot)


############################################


l = 2
h = 3
J = 1
n = 3*l*h
A_open = kagome(l,h; J=J, periodic=false); 
A      = kagome(l,h; J=J, periodic=true); 

val, X = QMC_SOC(A, verbose=false, pauli_NPA1=false)
plot = QMC_Plot(A_open,X)
name = "Kagome_" * string(l) * string(h) * "_" * string(J) * ".pdf"
Compose.draw(Compose.PDF(name, 16Compose.cm, 16Compose.cm), plot) #print to file
print(n); print("  "); print(val); plot

#######################################

# energy in VarBench normalization
nE = sum(A) / 2
println((val + nE/4)*4)

#######################################

println(nE)