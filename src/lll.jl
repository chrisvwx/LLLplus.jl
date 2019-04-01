"""
    B,T,Q,R = LLL(H,δ=3/4)

Do Lenstra–Lenstra–Lovász lattice reduction of matrix `H` using optional
parameter `δ`.  The output is `B`, an LLL-reduced basis; `T`, a unimodular
(meaning `det(T)=+/-1`) transformation matrix such that `B= H*T`; and
finally `Q` and `R` which are a QR decomposition of `B`.  So `H = B*inv(T) =
Q*R*inv(T)`.

Follows D. Wuebben, et al, "Lattice Reduction - A Survey with Applications
in Wireless Communications". IEEE Signal Processing Magazine, 2011. See
[`subsetsum`](@ref) for an application of `lll`.

# Examples
```jldoctest
julia> H= [1 2; 3 4];B,_ = lll(H); B
2×2 Array{Int64,2}:
 1  -1
 1   1

julia> H= BigFloat.([1.5 2; 3 4]) .+ 2im; B,_= lll(H); B
2×2 Array{Complex{BigFloat},2}:
 0.50+0.0im  0.0+1.0im
  1.0+0.0im  0.0+0.0im

julia> N=500;H = randn(N,N); B,T = lll(H);
```
"""
function lll(H::AbstractArray{Td,2},δ::Float64=3/4) where {Td<:Number}

if !(0.25 < δ < 1.0)
    error("δ must be between 1/4 and 1.");
end

B = copy(H);
L = size(B,2);
Qt,R = qr(B);
Q = Matrix(Qt); # A few cycles can be saved by skipping updates of the Q matrix.

# of course these should be replaced with different methods... I'll do it
# someday
if Td<:BigInt || Td<:BigFloat
    Ti = BigInt
elseif Td==Float32 || Td==Int32
    Ti=Int32
elseif Td==Int16
    Ti=Int16
elseif Td==Int8
    Ti=Int8
else
    Ti=Int
end
#Ti = (Td<:BigInt || Td<:BigFloat) ? BigInt : Int

if Td<:Complex
    T = Matrix{Complex{Ti}}(I, L, L)
else
    T = Matrix{Ti}(I, L, L)
end

lx  = 2;
while lx <= L

    # reduce lx-th column of B
    for k=lx-1:-1:1
        rk = R[k,lx]/R[k,k]
        mu = roundf(rk)
        if abs(mu)>zero(Ti)
            # B[:,lx]   -= mu * B[:,k]
            # R[1:k,lx] -= mu * R[1:k,k]
            # T[:,lx]   -= mu * T[:,k]
            B[:,lx]   .-= mu .* view(B,:,k)
            R[1:k,lx] .-= mu .* view(R,1:k,k)
            T[:,lx]   .-= mu .* view(T,:,k)
        end
    end

    nrm = norm(R[lx-1:lx,lx])
    if δ*abs(R[lx-1,lx-1])^2 > nrm^2

        # swap columns lx-1 and lx in B, T and R
        # B[:,[lx-1,lx]]    = B[:,[lx,lx-1]]
        # T[:,[lx-1,lx]]    = T[:,[lx,lx-1]]
        # R[1:lx,[lx-1,lx]] = R[1:lx,[lx,lx-1]]
        B[:,[lx-1,lx]]    .= view(B,:,[lx,lx-1])
        T[:,[lx-1,lx]]    .= view(T,:,[lx,lx-1])
        R[1:lx,[lx-1,lx]] .= view(R,1:lx,[lx,lx-1])

        # upper triangular by Givens rotation
        cc = R[lx-1,lx-1] / nrm # nrm = ||R[lx-1:lx,lx-1]|| after swapping
        ss = R[lx,lx-1]   / nrm
        Θ = [cc' ss; -ss cc]

        # the following multiply by Θ results in R[lx,lx-1] = 0
        R[lx-1:lx,lx-1:end] .= Θ * R[lx-1:lx,lx-1:end]
        # Q[:,lx-1:lx] = Q[:,lx-1:lx] * Θ'
        Q[:,lx-1:lx] .= view(Q,:,lx-1:lx) * Θ'
        lx = max(lx-1,2)
    else
        lx = lx+1;
    end
end
return B,T,Q,R
end


"""
    B = gauss(H)

Do Gauss (or Lagrange?!) reduction on the lattice defined by the two columns of
H.

Follows Fig 2.3 of "Lattice Basis Reduction: An Introduction to the LLL
Algorithm and Its Applications" by Murray R. Bremner, CRC Press, 2012.

# Examples
```jldoctest
julia> H = [1 2; 3 3]; B = gauss(H)
2×2 Array{Float64,2}:
 1.0  0.0
 0.0  3.0
```
"""
function gauss(H::AbstractArray{Td,2}) where Td
    if size(H,2)!=2
        error("Gauss reduction only works on two columns.")
    end
    x = float.(H[:,1])
    y = float.(H[:,2])
    if norm(x)<=norm(y)
        v1,v2=x,y
    else
        v1,v2=y,x
    end
    finished = false;
    while !finished
        m = round((v2'*v1)/(v1'*v1))
        v2 -= m*v1
        if norm(v1)<=norm(v2)
            finished = true
        else
            v1,v2 = v2,v1
        end
    end
    return [v1 v2]
end
