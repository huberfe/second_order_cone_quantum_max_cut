include("../QMC.jl");
include("../QMC_plot.jl");
include("../Lattice.jl");
include("../Lattice_other_square.jl")


println("Square 16 SOC")
A=Square(4, 4)
println("Objective")
t1=time()
println(QMC_SOC(A; pauli_NPA1=false, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")[1])
t2=time()
println("time elapsed:")
println(t2-t1)


println("Square 16 SOC+P1")
A=Square(4, 4)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))


println("Square 36 SOC+P1")
A=Square(6, 6)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))


println("Square 64 SOC+P1")
A=Square(8, 8)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))


println("Square 50 SOC+P1")
A=SquareRot()
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))

println("Square 100 SOC+P1")
A=Square(10, 10)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))


println("kagome 18 SOC+P1")
A=kagome(3, 2)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))


println("kagome 48 SOC+P1")
A=kagome(4, 4)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))



println("kagome 75 SOC+P1")
A=kagome(5, 5)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))



println("shuriken 24 SOC+P1")
A=shuriken(2, 2)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))



println("shuriken 96 SOC+P1")
A=shuriken(4, 4)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))



println("Pyrochlore 32 SOC+P1")
A=pyro_unit_4(2)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))



println("Pyrochlore 108 SOC+P1")
A=pyro_unit_4(3)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))


println("Pyrochlore 256 SOC+P1")
A=pyro_unit_4(4)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))



println("kagome 192 SOC+P1")
A=kagome(8, 8)
println("Objective")
val, X=QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=true, VarBench_norm=true, solver="scs")
println(val)
println("rounded objective")
println(rounded_obj(X, A, samps=1000, var_bench_norm=true))






