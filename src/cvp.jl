"""
    x=cvp(y,R,infinite=Val{true},Umax=1)

Solve the problem `argmin_x ||z-Hx||` for integer x using the technique from
the paper below, where H=QR and y=Q'*z. The input vector `y` is of length `n`,
with `H` of dimension `n` by `n`, and the returned vector `x` of length `n`. If
`infinite==Val{true}` then we search the (infinite) lattice, otherwise we
search integers in `[-Umax,Umax]`.  cvp does not handle complex
numbers.

Uses alg from "Faster Recursions in Sphere Decoding" Arash Ghasemmehdi, Erik
Agrell, IEEE Transactions on Information Theory, vol 57, issue 6 , June 2011.

# Examples
```jldoctest
julia> H=[1 2; 3 4]; Q,R=qr(H); uhat = cvp(Q'*[0,2],R)
2-element Array{Float64,1}:
  2.0
 -1.0

julia> n=100;H=randn(n,n);Q,R=qr(H);

julia> u=Int.(rand(0:1e10,n));y=H*u+rand(n)/100;

julia> uhat=cvp(Q'*y,R); sum(abs.(u-uhat))
0.0

julia> n=500;H=randn(n,n);Q,R=qr(H);

julia> u=Int.(rand(-1:1,n));y=H*u+rand(n)/10;

julia> uhat=cvp(Q'*y,R,Val{false},1); sum(abs.(u-uhat))
0.0
```
"""
function cvp(r::AbstractArray{Td,1},G::AbstractArray{Td,2},
               infinite=Val{true},Umax=1) where {Td<:Number}
    
    roundGA(x) = roundFinite(x,Td(Umax))
# nx=1
    n = length(r)
# mxNx=Int(ceil(log2(n)*100000))
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
    
    begin
        @label LOOP

#        nx+=1
        while λ[i]<C
            if i != 1
                i-=1
                for j=d[i]:-1:i+1
                    F[i,j-1] = F[i,j] - G[i,j]*u[j]
                end
                p[i] = F[i,i]/G[i,i]
                if infinite==Val{true}
                    u[i] = round(p[i])
                else
                    u[i] = roundGA(p[i])
                end
                y = (p[i]-u[i])*G[i,i]
                Δ[i] = signGA(y)
                λ[i] = λ[i+1]+y*y
            else
                û .= u
                C = λ[1]
            end
        end
        m = i

        while λ[i]≥C
#            if i==n || nx>mxNx
            if i==n
                # if nx>mxNx
                #     @warn("terminated search after $(nx) iterations.")
                # end
                return û
            else
                i+=1
                if infinite ≠ Val{true}
                    y = typemax(Td)
                end
                u[i]+=Δ[i]
                Δ[i]=-Δ[i]-signGA(Δ[i])
                if infinite ≠ Val{true}
                    if -Umax ≤ u[i] ≤ Umax
                        y = (p[i]-u[i])*G[i,i]
                    else
                        u[i]+=Δ[i]
                        Δ[i]=-Δ[i]-signGA(Δ[i])
                        if -Umax ≤ u[i] ≤ Umax
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

signGA(x) = x<=0 ? -one(x) : one(x)

function roundFinite(x,Umax)
    y=round(x)
    return  abs(y)>Umax ?  Umax*signGA(y) : y
end
    



"""
    b = svp(B)

Find the shortest basis vector `b` for the lattice formed by the matrix
`B`. This solves the 'shortest vector problem' (SVP). 

We use the [`cvp`](@ref) function in the library, which is not necessarily
the fastest SVP solver.  Roughly follows the CVP-to-SVP reduction in
http://web.eecs.umich.edu/~cpeikert/lic15/lec06.pdf

# Examples
```jldoctest
julia> H=[1 2; 3 4]; svp(H)
2-element Array{Int64,1}:
 -1
  1

julia> H= BigFloat.([2.5 2; 3 4]); svp(H)
2-element Array{BigFloat,1}:
  0.50
 -1.0 

```
"""
function svp(B::AbstractArray{Td,2}) where Td
    m,n= size(B)
    V = zeros(Td,m,n)
    c = zeros(Td,n)

    for i = 1:n
        Bi = copy(B)
        Bi[:,i] *=2
        Q,R = qr(Bi)
        vc = cvp(Q'*B[:,i],R)
        V[:,i] = Bi*vc-B[:,i]
        c[i] = V[:,i]'*V[:,i]
    end
    idx = sortperm(c)
    return V[:,idx[1]]
end

