include("../QMC.jl")
include("../QMC_plot.jl")
include("../Lattice.jl")
include("../pyrochlore_v2.jl")

A=pyro_unit_4(4)
mat_size=size(A)

println("trivial bound on objective")
println(-3*ones(Int64, 1, mat_size[1])*A*ones(Int64, mat_size[1], 1)/2)
###############################
t1=time()
println("basic one ")
println(QMC_SOC(A; pauli_NPA1=false, starbound=false, starbound4=false, SDP_form=false, verbose=false, VarBench_norm=true)[1])
print("time elapsed")
t2=time()
println(t2-t1)
############################
#println("now adding NPA1")
println(QMC_SOC(A; pauli_NPA1=true, starbound=false, starbound4=false, SDP_form=false, verbose=false, VarBench_norm=true)[1])
t3=time()
print("time elapsed")
t3-t2
