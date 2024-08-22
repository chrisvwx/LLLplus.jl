"""
    x=cvp(y,R::UpperTriangular)

Solve the problem `argmin_x ||y-Rx||` for integer x using the technique from
the paper below. The input vector `y` is of length `n`, with upper
triangular `R` of dimension `n` by `n`, and the returned vector `x` of
length `n`.

    x=cvp(y,R,infinite=Val(true),Umin=-1,Umax=-Umin,nxMax=Int(ceil(log2(n)*1e6)))

If `infinite==Val(true)` then we search the (infinite) lattice, otherwise we
search integers in `[Umin,Umax]`. When `nxMax` is included, it gives a
maximum number of steps to take, otherwise a heuristic value is used.
Note that `cvp` does not handle complex numbers.

Follows "Faster Recursions in Sphere Decoding" Arash Ghasemmehdi, Erik
Agrell, IEEE Transactions on Information Theory, vol 57, issue 6 , June 2011.

# Examples
```jldoctest
julia> H=[1 2; 3 4]; Q,R=qr(H); uhat = cvp(Q'*[0,2],UpperTriangular(R))
2-element Vector{Float64}:
  2.0
 -1.0

julia> uhat = cvp(Q'*[0,2],UpperTriangular(R),Val(false),0,100)
2-element Vector{Float64}:
 1.0
 0.0

julia> n=100;H=randn(n,n); Q,Rt=qr(H); R=UpperTriangular(Rt);

julia> u=Int.(rand(0:1e10,n));y=H*u+rand(n)/100;

julia> uhat=cvp(Q'*y,R); sum(abs.(u-uhat))
0.0

julia> n=500;H=randn(n,n); Q,Rt=qr(H); R=UpperTriangular(Rt);

julia> u=Int.(rand([-1,1],n));y=H*u+rand(n)/100;

julia> uhat=cvp(Q'*y,R,Val(false),-1,1); sum(abs.(u-uhat))
0.0
```
"""
function cvp(r::AbstractArray{Td,1},G::UpperTriangular{Td,Array{Td,2}},
             infinite=Val(true),Umin=-1,Umax=-Umin,
             nxMax=typemax(Int)) where {Td<:Number}

    roundGA(x) = roundFinite(x,Td(Umin),Td(Umax))
    n = length(r)
    C = Inf
    i = n+1
    d = ones(Int,n)*n
    λ = zeros(n+1)
    F = zeros(Td,n,n)
    F[:,n] = r
    p = zeros(n)
    u = NaN*zeros(Int,n)
    û = NaN*zeros(Int,n)
    Δ = zeros(Int,n)

    if nxMax == typemax(Int)
        nxMax=Int(ceil(log2(n)*1e6)) # heuristic
    end
    nx=1

    begin
        @label LOOP

        nx+=1
        while λ[i]<C
            if i != 1
                i-=1
                for j=d[i]:-1:i+1
                    F[i,j-1] = F[i,j] - G[i,j]*u[j]
                end
                p[i] = F[i,i]/G[i,i]
                if infinite==Val(true)
                    u[i] = round(p[i])
                else
                    u[i] = roundGA(p[i])
                end
                y = (p[i]-u[i])*G[i,i]
                Δ[i] = sgns(y)
                λ[i] = λ[i+1]+y*y
            else
                û .= u
                C = λ[1]
            end
        end
        m = i

        while λ[i]≥C
            if i==n || nx>nxMax
                if nx>nxMax
                    @warn("terminated search after $(nx) iterations.")
                end
                return û
            else
                i+=1
                if infinite ≠ Val(true)
                    y = typemax(Td)
                end
                u[i]+=Δ[i]
                Δ[i]=-Δ[i]-sgns(Δ[i])
                if infinite ≠ Val(true)
                    if Umin ≤ u[i] ≤ Umax
                        y = (p[i]-u[i])*G[i,i]
                    else
                        u[i]+=Δ[i]
                        Δ[i]=-Δ[i]-sgns(Δ[i])
                        if Umin ≤ u[i] ≤ Umax
                            y = (p[i]-u[i])*G[i,i]
                        end
                    end
                else
                    y = (p[i]-u[i])*G[i,i]
                end
                λ[i] = λ[i+1]+y*y
            end
        end
        d[m:i-1] .= i
        for j=m-1:-1:1
            if d[j]<i
                d[j]=i
            else
                @goto LOOP
            end
        end
        @goto LOOP
    end
end

sgns(x) = x<=0 ? -one(x) : one(x)

function roundFinite(x,Umin,Umax)
    y=round(x)
    if x<Umin
        return Umin
    elseif x>Umax
        return Umax
    else
        return y
    end
end



"""
    x=svp(B)

Find the shortest basis vector `b` for the lattice formed by the matrix
`B`. This solves the 'shortest vector problem' (SVP). Note that `svp`
does not handle complex numbers.

Follows "Algorithm SHORTESTVECTOR(G)" from "Closest Point Search in
Lattices" by Erik Agrell, Thomas Eriksson, Alexander Vardy, and
Kenneth Zeger in IEEE Transactions on Information Theory, vol. 48, no. 8,
August 2002.

This function IS BROKEN; see the second example below.

# Examples
```
julia> H=[1 2; 3 4]; LLLplus.svp(H)
2-element Vector{Int64}:
 -1
 -1

julia> H = [1 0 0 0;   0 1 0 0;   208 175 663 0;     651 479 0  663];

julia> LLLplus.svp(H)
4-element Vector{Int64}:
  16
 -19
   3
 -11

```
"""
function svp(G::AbstractArray{Td,2}) where {Td<:Number}

    @warn("This function appears broken; see the example...")

    m,n = size(G)
    if n==1 return G[:,1];    end

    G2,_ = lll(G)
    # if qr gave positive values on diagonal we'd use  Q,G3 = qr(G2)
    G3 = cholesky(Hermitian(G2'*G2)).U
    H3 = inv(G3)
    #Q = G2*H3
    uhat = decodeSVPAgrell(H3)
    yhat = G2*uhat # the "shortest vector"
    return yhat
end

"""
    z=svpu(B)

Find the vector of integers z such that B*z is the shortest vector
in the lattice with basis B.

# Examples
```
julia> H=[1 2; 3 4]; s=svp(H); u=LLLplus.svpu(H); [u s H*u]
2×3 Matrix{Int64}:
  1  -1  -1
 -1  -1  -1

```
"""
function svpu(G::AbstractArray{Td,2}) where {Td<:Number}
    m,n = size(G)
    if n==1 return [1]; end
    G2,T = lll(G)
    # if qr gave positive values on diagonal we'd use  Q,G3 = qr(G2)
    G3 = cholesky(Hermitian(G2'*G2)).U
    H3 = inv(G3)
    uhat = decodeSVPAgrell(H3)
    return T*uhat
end

"""
    x=decodeSVPAgrell(x,H)

Internal function.

Follows "Algorithm DECODE(G)" with "SHORTESTVECTOR" changes from the paper
"Closest Point Search in Lattices" by Erik Agrell, Thomas Eriksson,
Alexander Vardy, and Kenneth Zeger in IEEE Transactions on Information
Theory, vol. 48, no. 8, August 2002.


# Examples
```jldoctest
julia> H = inv([1.0 2; 0 4]);

julia> uhat = LLLplus.decodeSVPAgrell(H)
2-element Vector{Int64}:
 -1
  0


```
"""
function decodeSVPAgrell(H::AbstractArray{Td,2}) where {Td<:Number}

    n = size(H,1)
    bestdist=Inf
    k=n
    dist = zeros(n)
    e = zeros(n,n)
#    e[:,k] = H*x
    u = zeros(Int,n)
    û = zeros(Int,n)
    u[k] = round(e[k,k])
    y = (e[k,k]-u[k])/H[k,k]
    Δ = NaN*zeros(n)
    Δ[k] = sgns(y)
    begin
        @label LOOP
        newdist = dist[k]+y*y
        if newdist<bestdist
            if k!=1
                e[1:k-1,k-1]=e[1:k-1,k]-y*H[1:k-1,k]
                k -=1
                dist[k]=newdist
                u[k] = round(e[k,k])
                y = (e[k,k]-u[k])/H[k,k]
                Δ[k] = sgns(y)
            else
                if !(newdist≈0.0)
                    û[:] = u[:]
                    bestdist=newdist
                    k+=1
                end
                # û[:] = u[:]
                # bestdist=newdist
                # k+=1
                u[k] += Δ[k]
                y = (e[k,k]-u[k])/H[k,k]
                Δ[k] = -Δ[k] - sgns(Δ[k])
            end
        else
            if k==n
                return û
            else
                k+=1
                u[k]+= Δ[k]
                y = (e[k,k]-u[k])/H[k,k]
                Δ[k] = -Δ[k] - sgns(Δ[k])
            end
        end
        @goto LOOP
    end
end
