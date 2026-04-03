include("../QMC.jl");
include("../QMC_plot.jl");
include("../Lattice.jl");
include("../Lattice_other_square.jl")


println("Square 16 Perm Program")
A=Square(4, 4)
println("Objective")
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")


println("Square 16 Perm Program")
A=Square(4, 4)
println("Objective")
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")


println("Square 36 Perm Program")
A=Square(6, 6)
t1=time()
println("Objective")
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")


println("Square 64 Perm Program")
A=Square(8, 8)
println("Objective")
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")

println("Square 50 Perm Program")
A=SquareRot()
println("Objective")
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")

println("Square 100 SOC")
A=Square(10, 10)
println("Objective")
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")



println("kagome 18 Perm Program")
A=kagome(3, 2)
println("Objective")
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")





println("kagome 48 Perm Program ")
A=kagome(4, 4)
println("Objective")
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")


println("kagome 75 Perm Program ")
A=kagome(5, 5)
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")


println("shuriken 24 Perm Program")
A=shuriken(2, 2)
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")


println("shuriken 96 Perm Program")
A=shuriken(4, 4)
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")



println("Pyrochlore 32 Perm Program")
A=pyro_unit_4(2)
t1=time()
println(QMC_SOC(A; pauli_NPA1=false, starbound=false, starbound4=false, SDP_form=false, verbose=false, VarBench_norm=true, solver="scs")[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")


println("Pyrochlore 108 Perm")
A=pyro_unit_4(3)
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")


println("Pyrochlore 256 Perm")
A=pyro_unit_4(4)
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")



println("kagome 192 Perm")
A=kagome(8, 8)
t1=time()
println(l1_perm_program(A; solver="scs", verbose=false, VarBench_norm=false)[1])
t2=time()
println("time elapsed:")
println(t2-t1)
println("==============================================================")






