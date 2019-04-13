"""
    B, T = brun(H)

Brun's integer-relations alg implemented as a matrix decomposition. Takes as
input the matrix `H` and returns a reduced basis `B` and `T`, a unimodular
transformation matrix such that `B = H*T`. Brun reduction is often done with
`pinv(H)` as input to yield `B = pinv(H)*T`.

See V. Brun, "En generalisation av kjedebrøken I," Skr. Vid
ensk. Selsk. Kristiana, Mat. Nat. Klasse, 1919.  See
https://archive.org/stream/skrifterutgitavv201chri#page/300/mode/2up

Follows code from D. Wuebben, D. Seethaler, J. Jalden, and G. Matz, "Lattice
Reduction - A Survey with Applications in Wireless Communications" IEEE
Signal Processing Magazine, March 2011

# Examples
```jldoctest
julia> H=[1 2; 3 4]; B,T=brun(H); T
2×2 Array{Int64,2}:
  3  -1
 -2   1

```
"""
function brun(H::AbstractArray{Td,2}) where Td

    # get size
    K, M = size(H)

    if rank(H) < K
        error("Input basis does not have full row rank")
    end 

    Ti= getIntType(Td)    

    grammian = H'*H
    metric_vec = sqrt.(abs.(diag(grammian')))
    mu  = grammian[:,1] # initialize mu with an arbitrary column of the Grammian
    T = Matrix{Ti}(I, K, K)   # unimodular transformation matrix 
    B = copy(H) 

    doBrun = true
    while doBrun

        # Calculation of Brun's update value
        zw,s = findmax(mu)
        zw = copy(mu)
        zw[s] = -typemax(zw[1])
        zw,t = findmax(zw)
        r = round(mu[s]/mu[t])
        
        # basis update
        updated_vec = B[:,s] .- r'*B[:,t]
        
        # if metric can be improved
        updated_norm = norm(updated_vec)
        if updated_norm < metric_vec[s] 

            # update mu, unimodular and reduced matrices, metric
            mu[s] -= r*mu[t]
            T[:,s] .-= r'*view(T,:,t)  #T[:,s] -= r'*T[:,t]
            B[:,s] .= updated_vec
            metric_vec[s] = updated_norm

        else
            doBrun = false
        end

    end # 

    return B, T
end
