import Graphs
import GraphPlot
using Colors
import ColorSchemes
import NetworkLayout
using Compose

function QMC_Plot(A, X; large_edges=nothing)
    """ adjacency matrix A, edge values X, optional large_edges as list of (i,j) tuples """
    g = Graphs.SimpleGraph(A)

    X_max = maximum(X)
    Edges = collect(Graphs.edges(g))
    L = length(Edges)
    EdgeColors = zeros(RGB{Float64}, L)
    for it = 1:L
        i, j = Edges[it].src, Edges[it].dst
        EdgeColors[it] = get(
            reverse(ColorSchemes.colorschemes[:RdYlBu]),
            (X[i, j] + 1) / 2)
    end

    layout = (args...) -> GraphPlot.spring_layout(args...; C=1)

    # Compute node positions once
    locs_x, locs_y = layout(g)

    # Base plot with colored edges
    base_plot = GraphPlot.gplot(g,
        locs_x, locs_y,
        NODESIZE=0.02, nodefillc="black",
        edgestrokec=EdgeColors, EDGELINEWIDTH=1.5)

    if large_edges !== nothing && length(large_edges) > 0
        # Build a separate graph containing only the large edges
        n = Graphs.nv(g)
        g_overlay = Graphs.SimpleGraph(n)
        for (i, j) in large_edges
		if abs(A[i, j])>= 10^(-5)
	            Graphs.add_edge!(g_overlay, i, j)
		end
        end

        # Plot the overlay graph with black lines using the same positions
        L_overlay = Graphs.ne(g_overlay)
        overlay_plot = GraphPlot.gplot(g_overlay,
            locs_x, locs_y,
            NODESIZE=0.0,
            nodefillc=RGBA(0, 0, 0, 0),
            edgestrokec=[colorant"pink" for _ in 1:L_overlay],
            EDGELINEWIDTH=2.5,
            linetype="straight")

        # Apply dashing to the overlay by wrapping it with a Compose property
        dashed_overlay = compose(context(),
            (context(), Compose.strokedash([1mm, 1mm]),
             overlay_plot))

        # Compose overlay on top of base
        PLOT = compose(context(), dashed_overlay, base_plot)
    else
        PLOT = base_plot
    end

    return PLOT
end
