using Cairo, Compose, FileIO, JLD2
include("../Lattice.jl")
include("../QMC.jl")
include("../QMC_plot.jl")

function save_defects_data()
	# triangle graph with defects
	h = 12 #8  #use even numbers to avoid edges sticking out.
	l = 12 #12
	n = h*l

	g = TriangleGrid(l,h)
	A = Graphs.adjacency_matrix(g);

	#defect
	g_def = DeleteRandomEdges(g,0.99)  #0.97
	A_def = Graphs.adjacency_matrix(g_def);

	#GraphPlot.gplot(g, layout=GraphPlot.spectral_layout)
	#sprint layout (can be messed up too..)
	#layout=(args...)->GraphPlot.spring_layout(args...; C=1)
	#GraphPlot.gplot(g, layout=layout)

	# approximate max cut
	#val1, X1 = QMC_SOC(A,     pauli_NPA1=false, verbose=true);
	val2, X2 = QMC_SOC(A,     pauli_NPA1=true, verbose=true);
	obj, l_edges2=rounded_obj(X2, A)

	# with defects
	#val3, X3 = QMC_SOC(A_def, pauli_NPA1=false, verbose=true);
	val4, X4 = QMC_SOC(A_def, pauli_NPA1=true, verbose=true);
	obj, l_edges4=rounded_obj(X4, A)

	dict1=Dict("X2"=> X2, "X4"=> X4, "l_edges2"=>l_edges2, "l_edges4"=> l_edges4, "A"=>A, "g"=> g, "A_def"=> A_def, "g_def"=>g_def)
	print(l_edges2)
	print(l_edges4)

	FileIO.save("defect_vars"*string(round(Int, time()))*"point_99"*".jld2","dict1",dict1)
	FileIO.save("copy_defect_vars"*string(round(Int, time()))*"point_99"*".jld2","dict1",dict1)

	#plot1 = QMC_Plot(g,X1)
	#draw(PDF(join([n,"triangle_SOC.pdf"]), 16cm, 16cm), plot1)

	#plot2 = QMC_Plot(g,X2)
	#draw(PDF(join([n,"triangle_SOC_Pauli.pdf"]), 16cm, 16cm), plot2)

	#plot3 = QMC_Plot(g_def,X3)
	#draw(PDF(join([n,"triangle_def_SOC.pdf"]), 16cm, 16cm), plot3)

	#plot4 = QMC_Plot(g_def,X4)
	#draw(PDF(join([n,"triangle_def_SOC_Pauli.pdf"]), 16cm, 16cm), plot4)

end


function replot_defects()
	dict1=FileIO.load("defect_vars1771531439461584e9.jld2", "dict1")
	l=12
	h=12
	n=l*h

	g = TriangleGrid(l,h)
	A = Graphs.adjacency_matrix(g);


	X2=dict1["X2"]
	X4=dict1["X4"]
	l_edges2=dict1["l_edges2"]
	l_edges4=dict1["l_edges4"]
	A=dict1["A"]
	g=dict1["g"]
	A_def=dict1["A_def"]
	g_def=dict1["g_def"]




	samps=100

	for j=1:1:samps

		
		plot2 = QMC_Plot(g,X2, large_edges=l_edges2)
		draw(SVG(join([n,"triangle_SOC_Pauli"*string(j)*".svg"]), 16cm, 16cm), plot2)
	end
	for j=1:1:samps
		
		plot4 = QMC_Plot(g_def,X4, large_edges=l_edges4)
		draw(SVG(join([n,"triangle_def_SOC_Pauli"*string(j)*".svg"]), 16cm, 16cm), plot4)

	end




end
