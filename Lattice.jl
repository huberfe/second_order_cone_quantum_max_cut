# some periodic lattices:
# square
# shastry_sutherland2D
# ShastrySutherland2D_disorder
# ShastrySutherland2D_4RDM
# kagome
# shuriken
# TriangleGrid
# Erdos_renyi

import Graphs

function square(l,h;J=1, periodic=true)
    """ periodic 2D square lattice of length l and height h
        checked by hand on 3 x 3
    """
    @assert l > 2
    @assert h > 2
    
    L = l; H = h;
    n = L*H
    A = zeros(n,n)

    # indexing:
    # c = j*L + i
    # j counts rows up, i columns to the right
    # one row up: +L
    # one right:  +1
    
    # flattened cartesian coordinate, modulo torus.
    cc = (i,j) -> mod( mod(j+H, H)*L + mod(i+L,L), n) + 1

    for i in 0:L-1, j in 0:H-1        # i horizontal index, j vertical 
        # square grid: current coordinate
        c = cc(i,j)
        
        if i < L-1 || periodic
            A[c, cc(i+1,j)] = J   # horizontal edge right
        end
        if j < H-1 || periodic
            A[c, cc(i,j+1)] = J   # vertical edge up
        end
    end
    return A + A'
end

function shastry_sutherland2D(l,h; J=1, J_diag=1, periodic=true)
    """ periodic Shastry-Sutherland model in 2D,
        of  2*l (length/width) x 2*h (heigth)
        3x3 building blocks of the shape
        ---------
        |   | L |
        ---------
        | R |   |
        ---------
        R: diag right up / NE orientation
        L: diag left up / NW orientation

        periodic lattice only makes sense for l,h even!
        returns symmetric adjancency matrix with weights
    """

    # indexing:
    # c = j*L + i
    # j counts rows up, i columns to the right
    # one row up: +L
    # one right:  +1
    @assert l % 2 == 0
    @assert h % 2 == 0
    L = 2*l
    H = 2*h
    n = L*H
    A = zeros(n,n)

    # flattened cartesian coordinate, modulo torus.
    cc = (i,j) -> mod( mod(j+H, H)*L + mod(i+L,L), n) + 1

    for i in 0:L-1, j in 0:H-1        # i horizontal index, j vertical
        # square grid: current coordinate
        c = cc(i,j)

        if i < L-1 || periodic
            A[c, cc(i+1,j)] = J   # horizontal edge right
        end
        if j < H-1 || periodic
            A[c, cc(i,j+1)] = J   # vertical edge up
        end

        # only make diagonal-up links if not on last row
        if (j < H-1) || periodic
            # diag up right
            if mod(i,2) == 0 && mod(j,2) == 0 && (i<L-1 || periodic) # only if not on last column
                A[c, cc(i+1, j+1)] = J_diag # up right

            # diag up left
            elseif mod(i,2) == 0 && mod(j,2) == 1 && (0<i || periodic) #only if not in first column
                A[c, cc(i-1, j+1)] = J_diag # up left
            end
        end

    end
    #end if-else periodic
    return A + A'
end

function shastry_sutherland2D_4RDM(l,h)
    """ returns interesting subset of 4-RMDs
        (squares, horizontal/vertical nearest neighbor vertices)
        for periodic Shastry-Sutherland model in 2D,
        of  2*l (length/width) x 2*h (heigth)
        3x3 building blocks of the shape
        ---------
        |   | L |
        ---------
        | R |   |
        ---------
        R: diag right up / NE orientation
        L: diag left up / NW orientation

        periodic lattice only makes sense for l,h even!
        returns 4-RDMs along horizontal and vertical edges, and squares
    """
    # indexing:
    # c = j*L + i
    # j counts rows up, i columns to the right
    # one row up: +L
    # one right:  +1
    @assert l % 2 == 0
    @assert h % 2 == 0
    L = 2*l
    H = 2*h
    n = L*H
    #A = zeros(n,n)

    # flattened cartesian coordinate, modulo torus.
    cc = (i,j) -> mod( mod(j+H, H)*L + mod(i+L,L), n) + 1
    four_set = [] # set of 4-RDMs to conside: horizonal, vertical, square
    for i in 0:L-1, j in 0:H-1        # i horizontal index, j vertical
        # square grid: current coordinate
        c = cc(i,j)
        four_set = push!(four_set, sort([c, cc(i+1,j), cc(i+2,j), cc(i+3,j)]))    # horizonal line
        four_set = push!(four_set, sort([c, cc(i,j+1), cc(i,j+2), cc(i,j+3)]))    # vertical line
        four_set = push!(four_set, sort([c, cc(i,j+1), cc(i+1,j), cc(i+1,j+1)]))  # squares
    end
    return four_set
end


function shuriken(l,h; J=1, periodic=true)
    """ periodic Shuriken model in 2D,
        of  l (length/width) x h (heigth)
        1 building block of the shape
    """

    # indexing:
    # c = j*L + i
    # j counts rows up, i columns to the right
    # one row up: +L
    # one right:  +1
    # layout per unit cell
    #     4
    #  5 --- 3
    #  |     |  2
    #  0 --- 1

    L = l
    H = h
    U = 6 #elements in unit cell
    n = U*L*H
    A = zeros(n,n)

    # flattened cartesian coordinate, modulo torus.
    # i horizontal, j, k internal coordinate
    cc = (i,j,k=0) -> mod( U*( mod(j+H, H)*L + mod(i+L,L)) + k, n) + 1

    for i in 0:L-1, j in 0:H-1
        # square grid: current coordinate
        c = cc(i,j)

        #unit cell:
        #square
        A[c,c+1]   = J
        A[c+1,c+3] = J
        A[c+3,c+5] = J
        A[c+5,c]   = J

        #triangles to the right and up
        A[c+1,c+2]   = J
        A[c+2,c+3]   = J
        A[c+3,c+4]   = J
        A[c+4,c+5]   = J

        if i < L-1 || periodic
            A[c+2, cc(i+1,j)  ] = J   # cross right
            A[c+2, cc(i+1,j,5)] = J
        end
        if j < H-1 || periodic
            A[c+4, cc(i,j+1)  ]   = J   # cross up
            A[c+4, cc(i,j+1,1)]   = J
        end
    end

    return A + A'
end


function ShastrySutherland2D_disorder(l,h; J=1, J_diag=1, periodic=true, sig=0)
    """ periodic Shastry-Sutherland model in 2D, 
        of  2*l (length/width) x 2*h (heigth) 
        3x3 building blocks of the shape
        ---------
        |   | L |
        ---------
        | R |   |
        ---------
        R: diag right up / NE orientation
        L: diag left up / NW orientation

        periodic lattice only makes sense for l,h even!
        returns symmetric adjancency matrix with weights
    """

    # indexing:
    # c = j*L + i
    # j counts rows up, i columns to the right
    # one row up: +L
    # one right:  +1
    @assert l % 2 == 0
    @assert h % 2 == 0
    L = 2*l
    H = 2*h
    n = L*H
    A = zeros(n,n)

    # flattened cartesian coordinate, modulo torus.
    cc = (i,j) -> mod( mod(j+H, H)*L + mod(i+L,L), n) + 1

    for i in 0:L-1, j in 0:H-1        # i horizontal index, j vertical 
        # square grid: current coordinate
        c = cc(i,j)
        
        if i < L-1 || periodic
            A[c, cc(i+1,j)] = J*(1+sig*randn())   # horizontal edge right
        end
        if j < H-1 || periodic
            A[c, cc(i,j+1)] = J*(1+sig*randn())   # vertical edge up
        end

        # only make diagonal-up links if not on last row 
        if (j < H-1) || periodic
            # diag up right
            if mod(i,2) == 0 && mod(j,2) == 0 && (i<L-1 || periodic) # only if not on last column
                A[c, cc(i+1, j+1)] = J_diag*(1+sig*randn()) # up right

            # diag up left
            elseif mod(i,2) == 0 && mod(j,2) == 1 && (0<i || periodic) #only if not in first column
                A[c, cc(i-1, j+1)] = J_diag*(1+sig*randn()) # up left
            end
        end
        
    end
    #end if-else periodic
    return A + A'
end




function kagome(l,h; J=1, periodic=true)
    """ periodic kagome model in 2D,
        of  l (length/width) x h (heigth) 
        1 building block of the shape
    """

    # indexing:
    # c = j*L + i
    # j counts rows up, i columns to the right
    # one row up: +L
    # one right:  +1
    # layout per unit cell
    #    2
    #  /   \
    # 0 --- 1
    
    L = l
    H = h
    U = 3 #elements in unit cell
    n = U*L*H
    A = zeros(n,n)

    # flattened cartesian coordinate, modulo torus.
    # i horizontal, j, k internal coordinate
    cc = (i,j,k=0) -> mod( U*( mod(j+H, H)*L + mod(i+L,L)) + k, n) + 1 

    for i in 0:L-1, j in 0:H-1 
        # square grid: current coordinate
        c = cc(i,j)
        
        #unit cell:
        #triangle
        A[c,c+1]   = J
        A[c+1,c+2] = J
        A[c,c+2]   = J

        if i < L-1 || periodic
            A[c+1,cc(i+1,j,0)]  = J   #link right
        end
        
        if (j < H-1 || periodic)
            A[c+2,cc(i,j+1,0)]  = J   #link up (& right)
        end
        
        if (0 < i || periodic) && (j < H-1 || periodic)
            A[c+2,cc(i-1,j+1,1)] = J
        end
    end

    return A + A'
end


##############
# triangle lattices, with deformations

function TriangleGrid(l,h)
    """ triangle grid of side lenght l, height h
    """
    n = h*l
    g = Graphs.SimpleGraph(n)
    for i in 1:l, j in 1:h # i horizontal index, j vertical
        c = (j-1)*l+i  # current coordinate
        if  i != l # avoid right border
            Graphs.add_edge!(g, c, c+1) # horizontal edge
        end
        if j != h #avoid top border
            Graphs.add_edge!(g, c, c+l) # vertical edge edge
        end
        if j != h &&  i != l && mod(j,2) == 1 
            Graphs.add_edge!(g, c, c+l+1) # vertical edge edge # diag right up
        end
        if j != h &&  i != 1 && mod(j,2) == 0 
            Graphs.add_edge!(g, c, c+l-1) # vertical edge edge # diag left up
        end
    end    
    return g
end

function DeleteRandomEdges(g,p)
    """ delete random edges of the graph g with probability p.
        Careful, might produce disconnected graph!
    """
    h = copy(g)
    for e in Graphs.edges(h)
        if rand() > p
            Graphs.rem_edge!(h, e)
        end
    end
    return h
end

function erdos_renyi(n, p)
    #p is probability of edge formation
    A=zeros(n, n)
    for i=1:n, j=i+1:n
        if rand()<= p
            A[i, j]=1
        end
    end
    return A+A'
end

##helper scripts for making pyrochlore lattice:
##########################################
function ad_l(a, b, m)
    if m==0
        return a+b
    else
        return (a+b) .% m
    end
end

function get_ind(i, j, k, l, lx)
    return i*4*lx^2+4*lx*j+4*k+l+1
end
##########################################

function pyro_unit_4(lx)
    ##produces perdiodic pyrochlore lattice with unit cells of size 4 of size lx by lx by lx
    
    g = Graphs.SimpleGraph((lx-1)*4*lx^2+4*lx*(lx-1)+4*(lx-1)+3+1)

    vert_list=[]
    ind_list=[]
    for i in collect(0:1:lx-1)
        for j in collect(0:1:lx-1)
            for k in collect(0:1:lx-1)
                for l in collect(0:1:3)

                    push!(vert_list, [[i, j, k], l])
                    push!(ind_list, i*4*lx^2+4*lx*j+4*k+l+1)
                end
            end
        end
    end

    #return ind_list
    for i in collect(0:1:lx-1)
        for j in collect(0:1:lx-1)
            for k in collect(0:1:lx-1)
            
                cur_big_vert=[i, j, k]
                cur_ind0=get_ind(i, j, k, 0, lx)
                cur_ind1=get_ind(i, j, k, 1, lx)
                cur_ind2=get_ind(i, j, k, 2, lx)
                cur_ind3=get_ind(i, j, k, 3, lx)

                

                #plx_ind = vert_list.index([ad_l(cur_big_vert, [1, 0, 0], m=lx), 0])
                #ad_mat[cur_ind1, plx_ind]=1
                #ad_mat[plx_ind, cur_ind1]=1

                plx_vec=ad_l(cur_big_vert, [1, 0, 0], lx)
                plx_ind = get_ind(plx_vec[1], plx_vec[2], plx_vec[3], 0, lx)
                Graphs.add_edge!(g, cur_ind1, plx_ind)
                
                #ply_ind = vert_list.index([ad_l(cur_big_vert, [0, 1, 0], m=lx), 0])
                #ad_mat[cur_ind2, ply_ind]=1
                #ad_mat[ply_ind, cur_ind2]=1

                ply_vec=ad_l(cur_big_vert, [0, 1, 0], lx)
                ply_ind=get_ind(ply_vec[1], ply_vec[2], ply_vec[3], 0, lx)
                Graphs.add_edge!(g, cur_ind2, ply_ind)
        
                
                #plz_ind=vert_list.index([ad_l(cur_big_vert, [0, 0, 1], m=lx), 0])
                #ad_mat[cur_ind3, plz_ind]=1
                #ad_mat[plz_ind, cur_ind3]=1

                plz_vec=ad_l(cur_big_vert, [0, 0, 1], lx)
                plz_ind=get_ind(plz_vec[1], plz_vec[2], plz_vec[3], 0, lx)
                Graphs.add_edge!(g, cur_ind3, plz_ind)



                ##now need to add internal edges
                internal_list=[cur_ind0, cur_ind1, cur_ind2, cur_ind3]
                for q1 in internal_list
                    for q2 in internal_list
                        if q1!= q2
                            Graphs.add_edge!(g, q1, q2)
                        end
                    end
                end

                ##now other ones
                
                #ind=vert_list.index([ad_l(cur_big_vert, [0, -1, 1], m=lx), 3])
                #ad_mat[cur_ind2, ind]=1
                #ad_mat[ind, cur_ind2]=1

                vec=ad_l(cur_big_vert, [0, -1, 1], lx)
                ind=get_ind(vec[1], vec[2], vec[3], 3, lx)
                Graphs.add_edge!(g, cur_ind2, ind)

                #ind=vert_list.index([ad_l(cur_big_vert, [1, -1, 0], m=lx), 3])
                #ad_mat[cur_ind1, ind]=1
                #ad_mat[ind, cur_ind1]=1
                vec=ad_l(cur_big_vert, [1, -1, 0], lx)
                ind=get_ind(vec[1], vec[2], vec[3], 3, lx)
                Graphs.add_edge!(g, cur_ind1, ind)

                #ind=vert_list.index([ad_l(cur_big_vert, [1, 0, -1], m=lx), 2])
                #ad_mat[cur_ind1, ind]=1
                #ad_mat[ind, cur_ind1]=1
                vec=ad_l(cur_big_vert, [1, 0, -1], lx)
                ind=get_ind(vec[1], vec[2], vec[3], 2, lx)
                Graphs.add_edge!(g, ind, cur_ind1)


            end
        end
    end

    return g
end

