#since we only are comparing this to a single square lattice Ill write it for the specific size

#there are different parameterizations oif the points in the grid
#index is the order which they will show up in the matrix
#(k, l) corresponds to the position in the grid where k is the row and l is the column.  k and l range from 0 to 9 but not every
#point in {0, ..., 9} x {0, ..., 9} is a lattice point.  In the first column the points go as (0, 1), (0, 3), (0, 5), (0, 7) and (0, 9).  
#In the second column it is (1, 0), (1, 2), (1, 4), (1, 6), (1, 8).  In the third is is (2, 1), (2, 3), (2, 5), (2, 7), and (2, 9) etc.
#In the grid a point (k, l) is connected to (k, l)+(1, 1), (k, l)+(-1, -1), (k, l)+(1, -1), (k, l)+(-1, 1).

#theres also an intermediate grid Im using to go from index to grid.  In this we simply have 10 rows and 5 columns with the points indexed 
#in the obvious way.  In this grid (i, j) is mapped to index 5i+j+1


n=50

function grid_to_ind(k, l)
	i=k
	if mod(k, 2)==1
		j=Int(round(l/2))
	end
	if mod(k, 2)==0
		
		j=Int(round((l-1)/2))
	end
	return 5*i+j+1
end

function ind_to_grid(ind)
	j=mod(ind-1, 5)
	i=Int(round((ind-j-1)/5))
	
	k=i 
	if mod(i, 2)==1
		l=2*j
	end
	if mod(i, 2)==0
		l=2*j+1
	end
	return (k, l) 
end



function SquareRot()
	##points will be indexed as (i, i+1) or (i+1, i) mod h
	#(i, j) will correspond to index 10
	A=zeros(50, 50)
	for i=1:50
		temp=ind_to_grid(i)
		k=temp[1]
		l=temp[2]

		#println(i)
		#println(grid_to_ind(k, l))
		
		k1=mod(k+1, 10)
		l1=mod(l+1, 10)

		k2=mod(k+1, 10)
		l2=mod(l-1, 10)

		#println("==============================================")
		#println(temp)
		#println(k1)
		#println(l1)
		#println(k2)
		#println(l2)

		ind1=grid_to_ind(k1, l1)
		ind2=grid_to_ind(k2, l2)

		A[i, ind1]=1
		A[i, ind2]=1

	end
	return A+A'

end

#println(ind_to_grid(2))
#println(grid_to_ind(0, 3))

#println(SquareRot()*ones(50, 1))