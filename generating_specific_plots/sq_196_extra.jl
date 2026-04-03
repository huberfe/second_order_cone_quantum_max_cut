include("../QMC.jl");
include("../QMC_plot.jl");
include("../Lattice.jl");
include("../Lattice_other_square.jl")

println("Square 16 SOC+P1")
A=Square(4, 4)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))
println("============================================")

println("Square 196 SOC+P1")
A=Square(14, 14)
println("Objective")
t1=time()
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
t2=time()
println("time elapsed:")
println(t2-t1)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))
