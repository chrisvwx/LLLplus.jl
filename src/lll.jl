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

Ti= getIntType(Td)
T = Matrix{Ti}(I, L, L)
zeroTi = real(zero(Ti))

lx  = 2;
while lx <= L

    # reduce lx-th column of B
    for k=lx-1:-1:1
        rk = R[k,lx]/R[k,k]
        mu = roundf(rk)
        if abs(mu)>zeroTi
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

Do Gauss/Lagrange reduction on the lattice defined by the two columns of
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


"""
    getIntType(Td)

    Return an integer type that matches Td.  I.e. BigInt for Td=BigFloat, Int64
    for Td=Float64. Should replace this with multiple functions; the
    following gives a hint, but it's not correct:

getIntType(Td::Complex{Tr}) = Complex{getIntType(Tr)}
getIntType(Td) where {Td<:Integer} = Td
getIntType(<:Float64) = Int64
getIntType(<:Float32) = Int32
getIntType(<:Float16) = Int16
getIntType(<:BigFloat) = BigInt

    this is better:
getIntType(Td::Type{Tr}) where {Tr<:Float128} = Int128

"""
function getIntType(Td)
    Tr = real(Td)
    if Tr<:BigInt || Tr<:BigFloat
        Ti = BigInt
    elseif Tr==Float32 || Tr==Int32
        Ti=Int32
    elseif Tr==Int16
        Ti=Int16
    elseif Tr==Int8
        Ti=Int8
    else
        Ti=Int
    end
    if Td<:Complex
        Ti = Complex{Ti}
    end
    return Ti
end

"""
    B,T,Q,R = sizereduction(H)
Do size reduction of matrix H.  The output is B, an size-reduced
basis; T, a unimodular (det(T)=+/-1) transformation matrix such that B=
H*T; Q and R such that B=Q*R and Q is orthonormal, and R is upper
triangular.  So H = Q*R*inv(T)

H can be of Integer, FloatingPoint, BigInt, or BigFloat types. The core
algorithm is designed for floating-point.

Size reduction is the first part of LLL lattice reduction.

# Examples
```jldoctest
julia> H= [1 2; 3 4];B,_ = sizereduction(H); B
2×2 Array{Int64,2}:
 3   1
 1   1

julia> H= BigFloat.([1.5 2; 3 4]) .+ 2im; B,_= sizereduction(H); B
2×2 Array{Complex{BigFloat},2}:
 1.5+2.0im  0.50+0.0im
 3.0+2.0im   1.0+0.0im

julia> N=100;H = randn(N,N); B,T = sizereduction(H);

```
"""
function sizereduction(H::Array{Td,2}) where {Td}
    B = copy(H);
    L = size(B,2);
    Q,R = qr(B);

    Ti= getIntType(Td)
    T = Matrix{Ti}(I, L, L)

    for lx=2:L
        # reduce lx-th column of B
        for k=lx-1:-1:1
            rk = R[k,lx]/R[k,k]
            mu = roundf(rk)
            if abs(mu)>0
                # B[:,lx]   -= mu * B[:,k]
                # R[1:k,lx] -= mu * R[1:k,k]
                # T[:,lx]   -= mu * T[:,k]
                B[:,lx]   .-= mu .* view(B,:,k)
                R[1:k,lx] .-= mu .* view(R,1:k,k)
                T[:,lx]   .-= mu .* view(T,:,k)
            end
        end
    end

    return B,T,Q,R
end
