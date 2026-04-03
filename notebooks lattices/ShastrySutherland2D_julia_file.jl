include("../QMC.jl")
include("../QMC_plot.jl")
include("../Lattice.jl")

import Cairo # Cairo is needed for Compose.PDF
import Compose #import Cairo

##########################################

get(reverse(ColorSchemes.colorschemes[:RdYlBu]), 0)
reverse(ColorSchemes.colorschemes[:diverging_linear_bjy_30_90_c45_n256])
reverse(ColorSchemes.colorschemes[:RdYlBu])

##########################################


#Check correctness of lattice
A = ShastrySutherland2D(2,2;J=1,J_diag=-1, periodic=false)
g = Graphs.SimpleGraph(A)
import GraphPlot
plot=GraphPlot.gplot(g, layout=GraphPlot.spectral_layout)
Compose.draw(Compose.PDF("SS_model_fig1.pdf", 16Compose.cm, 16Compose.cm), plot) #print to file

##########################################

l = h = 4
J_diag = 1
# open boundary lattice only for plotting
A_open = ShastrySutherland2D(l,h; J=1, J_diag=1, periodic=false); 

##########################################

J = 0.45

print(J," / ", J_diag, "  ")
A = ShastrySutherland2D(l,h; J=J, J_diag=J_diag, periodic=true);
val, X = QMC_SOC(A, verbose=false, pauli_NPA1=false)
plot = QMC_Plot(A_open,X)
name = "ShastrySutherland_" * string(l) * string(h) * "_" * string(J) * ".pdf"
Compose.draw(Compose.PDF(name, 16Compose.cm, 16Compose.cm), plot) #print to file
print(val); #plot

##########################################


#test: same but with selected 4-RDMs (too large for my laptop)

J = 0.45

FS = ShastrySutherland2D_4RDM(l,h)
print(J," / ", J_diag, "  ")
A = ShastrySutherland2D(l,h; J=J, J_diag=J_diag, periodic=true);
val, X = QMC_SOC(A, verbose=false, pauli_NPA1=false, four_qb = FS)
plot = QMC_Plot(A_open,X)
name = "ShastrySutherland_" * string(l) * string(h) * "_" * string(J) * "4RDM" * ".pdf"
Compose.draw(Compose.PDF(name, 16Compose.cm, 16Compose.cm), plot) #print to file
print(val); #plot

##########################################


J = 0.5

print(J," / ", J_diag, "  ")
A = ShastrySutherland2D(l,h; J=J, J_diag=J_diag, periodic=true);
val, X = QMC_SOC(A, verbose=false, pauli_NPA1=false)
plot = QMC_Plot(A_open,X)
name = "ShastrySutherland_" * string(l) * string(h) * "_" * string(J) * ".pdf"
Compose.draw(Compose.PDF(name, 16Compose.cm, 16Compose.cm), plot) #print to file
print(val); plot

##########################################

J = 0.6

print(J," / ", J_diag, "  ")
A = ShastrySutherland2D(l,h; J=J, J_diag=J_diag, periodic=true);
val, X = QMC_SOC(A, verbose=false, pauli_NPA1=false)
plot = QMC_Plot(A_open,X)
name = "ShastrySutherland_" * string(l) * string(h) * "_" * string(J) * ".pdf"
Compose.draw(Compose.PDF(name, 16Compose.cm, 16Compose.cm), plot) #print to file
print(val); #plot

##########################################

# phase transition dimer - plaquette: J/J_D = 0.72

l = 4
h = 4

# here: exactly at 0.5 the plaquette-Neel transition happens
J = 0.5
J_diag = 1

A      = ShastrySutherland2D(l,h; J=J, J_diag=J_diag, periodic=true)
A_open = ShastrySutherland2D(l,h; J=J, J_diag=J_diag, periodic=false)

val, X = QMC_SOC(A; pauli_NPA1=true, verbose=false)
print(J," / ", J_diag)
plot = QMC_Plot(A_open, X)
name = "ShastrySutherland_" * string(l) * string(h) * "_" * string(J) * "PauliNPA" * ".pdf"
Compose.draw(Compose.PDF(name, 16Compose.cm, 16Compose.cm), plot) #print to file
#plot

##########################################

# phase transition dimer - plaquette: J/J_D = 0.72

l = h = 5

# here: exactly at 0.5 the plaquette-Neel transition happens
J = 0.51
J_diag = 1

A      = ShastrySutherland2D(l,h; J=J, J_diag=J_diag, periodic=true)
A_open = ShastrySutherland2D(l,h; J=J, J_diag=J_diag, periodic=false)

val, X = QMC_SOC(A; pauli_NPA1=true, verbose=false)
print(J," / ", J_diag)
plot = QMC_Plot(A_open, X)
name = "ShastrySutherland_" * string(l) * string(h) * "_" * string(J) * "PauliNPA" * ".pdf"
draw(Compose.PDF(name, 16cm, 16cm), plot)
#plot