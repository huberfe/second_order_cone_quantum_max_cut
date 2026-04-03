include("../QMC.jl");
include("../QMC_plot.jl");
include("../Lattice.jl");
include("../Lattice_other_square.jl")


###first the funny rotated square
A=SquareRot()
println("first thing is the square 50")
t1=time()
println("basic one ")
println(QMC_SOC(A; pauli_NPA1=false, starbound=false, starbound4=false, SDP_form=false, verbose=false, VarBench_norm=true)[1])
print("time elapsed")
t2=time()
println(t2-t1)
############################
println("now adding NPA1")
println(QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=false, VarBench_norm=true)[1])
t3=time()
print("time elapsed")
t3-t2

println("===================================")
println("now the basic realxation for square 196")
A=Square(14,14;J=1, periodic=true)
t1=time()
println("basic one ")
println(QMC_SOC(A; pauli_NPA1=false, starbound=false, starbound4=false, SDP_form=false, verbose=false, VarBench_norm=true)[1])
print("time elapsed")
t2=time()
println(t2-t1)

println("===================================")
println("Now basic relaxation for Kagome 192")
A=kagome(8,8; J=1, periodic=true)
t1=time()
println("basic one ")
println(QMC_SOC(A; pauli_NPA1=false, starbound=false, starbound4=false, SDP_form=false, verbose=false, VarBench_norm=true)[1])
print("time elapsed")
t2=time()
println(t2-t1)

println("===================================")
println("now the advanced realxation for square 196")
A=Square(14,14;J=1, periodic=true)
t1=time()
println("basic one ")
println(QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=false, VarBench_norm=true)[1])
print("time elapsed")
t2=time()
println(t2-t1)

println("===================================")
println("Now advanced relaxation for Kagome 192")
A=kagome(8,8; J=1, periodic=true)
t1=time()
println("basic one ")
println(QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=false, VarBench_norm=true)[1])
print("time elapsed")
t2=time()
println(t2-t1)
