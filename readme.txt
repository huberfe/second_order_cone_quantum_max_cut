Code for:

Second order cone relaxations for quantum Max Cut
Felix Huber, Kevin Thompson, Ojas Parekh, Sevag Gharibian
https://arxiv.org/abs/2411.04120




Main files:

    QMC.jl
    QMC_plot.jl
    Lattice.jl
    SwapMat.jl

Requirements:

    JuMP, LinearAlgebra, SparseArrays
    Graphs, GraphPlot, Colors, ColorSchemes, NetworkLayout, Compose


###############################

Main functions:

QMC.jl: 

QMC_SOC(A; pauli_NPA1=false, four_qb=false, partition_31=true, partition_31_relax=false,
                    starbound=false, starbound4=false, SDP_form=false,
                    solver="mosek", verbose=false, VarBench_norm=false)
    returns both lower bound and edge values from SOC/SDP

HamiltonianMaxCut(A; is_sp=false::Bool)
    returns exact ground state energy by diagonalization

QMC_to_varbench(A, val)
    """ transforms QMC to Variational Benchmark normalization """

rounded_obj(X, A; samps=1000, var_bench_norm=true)
returns a feasible solution after rounding + list of large edges.




QMC_Plot.jl:

QMC_Plot(A, X; large_edges=nothing)
    """ adjacency matrix A, edge values X, optional large_edges as list of (i,j) tuples """
    plots graph with edge values, optionally with large values overlaid.
