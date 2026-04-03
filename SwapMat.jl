vec2num(v, b) = foldl((x,y) -> b*x+y, v; init=0)

function swap_padded(p::Int, q::Int, n::Int; d=2::Int, is_sp=false::Bool)
    if is_sp
        res = SparseArrays.spzeros(d^n, d^n)
        v = zeros(Int, n)
        for ind in 0:d^n-1
            digits!(v, ind; base=d)
            reverse!(v)
            v[p], v[q] = v[q], v[p]
            res[ind+1, vec2num(v, d)+1] = 1 
        end
        return res

    else
        """ code from Sebastian Designolle"""
        res = zeros(d^n, d^n)
        v = zeros(Int, n)
        for ind in 0:d^n-1
            digits!(v, ind; base=d)
            reverse!(v)
            v[p], v[q] = v[q], v[p]
            res[ind+1, vec2num(v, d)+1] = 1 
        end
        return res
    end


#    """ code from Sebastian Designolle"""
 #   res = zeros(d^n, d^n)
  #  v = zeros(Int, n)
   # for ind in 0:d^n-1
    #    digits!(v, ind; base=d)
     #   reverse!(v)
      #  v[p], v[q] = v[q], v[p]
       # res[ind+1, vec2num(v, d)+1] = 1 
    #end
    #return res
end

function swap(d=2)
    S = zeros(d^2,d^2)
    for i=1:d, j=1:d  # |ij><ji|
        S[i+d*(j-1),j+d*(i-1)] = 1
    end
    return S
end

@assert swap_padded(1,2,2) == swap(2)
