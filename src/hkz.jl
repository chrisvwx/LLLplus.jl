"""
    B,T=hkz(H)

Do Hermite-Korkine-Zolotarev (HKZ) reduction of the basis `H`, returning
the reduced basis `B` and the unimodular rotation `T` such that B = H*T.
HKZ reduction is sometimes called "Hermite-Korkine-Zolotareff" or
"Korkine-Zolotareff" reduction.

Based on "Practical HKZ and Minkowski Lattice Reduction Algorithms" by Wen
Zhang, Sanzheng Qiao, and Yimin Wei, 17 Aug. 2011.
http://www.cas.mcmaster.ca/~qiao/publications/ZQW11.pdf

# Examples
```jldoctest
julia> H=[1 2; 3 4];

julia> B,T= LLLplus.hkz(H); B
2×2 Matrix{Int64}:
  1  -1
 -1  -1

julia> N=9; H=rand(0:10,N,N); B,_=LLLplus.hkz(H);

julia> A = [10.6347 -66.2715  9.3046 17.5349 24.9625 # From
                0     8.6759 -4.7536 -3.9379 -2.3318 # TransIT, v65,
                0        0    0.3876  0.1296 -0.2879 # n3, Mar 2019,
                0        0    0       0.0133 -0.0082 # p 1929.
                0        0    0       0       0.0015];

julia> B,_= LLLplus.hkz(A); # our HKZ doesn't work well on A ...

```
"""
function hkz(B::AbstractArray{Td,2}) where {Td<:Number}
    Z,_ = hkz_red(B)
    return B*Z,Z
end


"""
    x=hkz_red(B)

Do Hermite-Korkine-Zolotareff reduction of the basis B.

This is the "HKZ-red" function for HKZ reduction in
"Practical HKZ and Minkowski Lattice Reduction Algorithms" by
Wen Zhang, Sanzheng Qiao, and Yimin Wei, 17 Aug. 2011.
http://www.cas.mcmaster.ca/~qiao/publications/ZQW11.pdf

# Examples
```jldoctest
julia> B=[1 2; 3 4];

julia> Z,_= LLLplus.hkz_red(B);

julia> B*Z
2×2 Matrix{Int64}:
  1  -1
 -1  -1

julia> N=10; Bo=rand(0:100,N,N); Z,_=LLLplus.hkz_red(Bo); B=Bo*Z;

julia> A = [10.6347 -66.2715  9.3046 17.5349 24.9625 # From
                0     8.6759 -4.7536 -3.9379 -2.3318 # TransIT, v65,
                0        0    0.3876  0.1296 -0.2879 # n3, Mar 2019,
                0        0    0       0.0133 -0.0082 # p 1929.
                0        0    0       0       0.0015];

julia> Z,_= LLLplus.hkz_red(A); # our HKZ doesn't work well on A ...

```
"""
function hkz_red(B::AbstractArray{Td,2}) where {Td<:Number}
    m,n = size(B)
    Ti= getIntType(Td)

    #Q,R = qr(B)  # QR doesn't give positive diagonal elements
    # HKZ-Red doesn't specify positive diagonal, but elsewhere in the
    # paper it does for QR, so I'll use it
    R = Matrix(cholesky(B'*B).U)
    Z = Matrix{Ti}(I, n,n)
    for k=1:n-1
        z = svpu(R[k:n,k:n])
        transform!(R,Z,z,k)
    end
    # size reduce
    for lx=2:n
        # reduce lx-th column of B
        for k=lx-1:-1:1
            rk = R[k,lx]/R[k,k]
            mu = round(rk)
            if abs(mu)>0
                R[1:k,lx] .-= mu .* view(R,1:k,k)
                Z[:,lx]   .-= mu .* view(Z,:,k)
            end
        end
    end
    return Z,R
end


"""
    transform!(R,Z,z,k)

Basis expansion. Given n x n upper triangular R, n x n unimodular Z, n-k+1
vector of integer z, and index k, return updated R and Z such that
R[k,k] = ||R[k:n,k:n]*z||.

This is the "TRANSFORM" basis expansion function from
"Practical HKZ and Minkowski Lattice Reduction Algorithms" by
Wen Zhang, Sanzheng Qiao, and Yimin Wei, 17 Aug. 2011.
http://www.cas.mcmaster.ca/~qiao/publications/ZQW11.pdf

"""
function transform!(R::AbstractArray{Td,2},Z::AbstractArray{<:Integer,2},
                    z::AbstractArray{<:Integer,1},k::Integer) where {Td<:Number}
    n = size(R,1)
    for j=n-k+1:-1:2
        if z[j]≠0
            # d=gcd(z[j-1],z[j]),  a,b are such that a*z[j−1] + b*z[j] = d
            d,a,b=gcdx(z[j-1],z[j])
            M = [z[j-1]÷d -b;
                 z[j]÷d    a]
            z[j-1] = d
            R[1:j+k-1, j+k-2:j+k-1]= R[1:j+k-1, j+k-2:j+k-1]*M
            c,s = givensrotation(R[j+k-2,j+k-2],R[j+k-1,j+k-2])
            G = [c -s; s c]
            R[j+k-2:j+k-1,j+k-2:n] = G'*R[j+k-2:j+k-1,j+k-2:n]
            Z[:,j+k-2:j+k-1] = Z[:,j+k-2:j+k-1]*M
        end
    end
end


"""
    c,s = givensrotation(a,b)
"""
function givensrotation(a,b)
    if b == 0
        c = 1
        s = 0
    else
        if abs(b) > abs(a)
            r = a / b
            s = 1 / sqrt(1 + r^2)
            c = s*r
        else
            r = b / a
            c = 1 / sqrt(1 + r^2)
            s = c*r
        end
    end
    return c,s
end


"""
    ishkzreduced(B)

Determine if the matrix B is Hermite-Korkine-Zolotarev (HKZ) reduced or
not. See the `hkz` function.

# Examples
```
julia> H= [1 2; 3 4];ishkzreduced(H)
false

julia> B,_=hkz(H); ishkzreduced(B)
true

```
"""
function ishkzreduced(B)
    if issizereduced(B)
        m,n = size(B)
        R = Matrix(cholesky(B'*B).U)
        for i=1:n
            s = svp(R[i:n,i:n])
            if !(R[i,i]≈sqrt(s'*s))
                return false
            end
        end
        return true
    else
        return false
    end
end
