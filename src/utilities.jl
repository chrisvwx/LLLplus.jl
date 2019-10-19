"""
    orthogonalitydefect(B)

Find the orthogonality defect of the matrix B defined, for example,
on page 2 of Bennet 

[Bennet](http://users.eecs.northwestern.edu/~hbennett/publications/bases.pdf)

# Examples
```jldoctest
julia> H= [1 2; 3 4];B,T=lll(H);

julia> [orthogonalitydefect(H) orthogonalitydefect(B)]
1×2 Array{Float64,2}:
 7.07107  1.0

```
"""
function orthogonalitydefect(B)
    N = size(B,2)
    prod = one(B[1,1])
    for n=1:N
        prod *= sqrt(B[:,n]'*B[:,n])
    end
    return prod/abs(det(B))
end




"""
    hermitefactor(B)

Find the Hermite factor of matrix B

# Examples
```jldoctest
julia> H= [1 2; 3 4];hermitefactor(H)
1.5811388300841898

```
"""
hermitefactor(B) = sqrt(B[:,1]'*B[:,1])/abs(det(B))




"""
    seysencond(B)

Seysen condition number as on, for example, page 3 of Bennet 

[Bennet](http://users.eecs.northwestern.edu/~hbennett/publications/bases.pdf)

# Examples
```jldoctest
julia> H= [1 2; 3 4];seysencond(H)
2.8284271247461903

```
"""
function seysencond(B)
    N = size(B,2)
    Bi = inv(B)'
    quotients = zeros(typeof(Bi[1,1]),N)
    for n=1:N
        quotients[n] = sqrt(B[:,n]'*B[:,n])/sqrt(Bi[:,n]'*Bi[:,n])
    end
    return maximum(quotients)
end



"""
    islllreduced(B)

Determine if the matrix B is LLL reduced or not. See p 56 of Bremner for a
definition. 

M. R. Bremner, "Lattice Basis Reduction: An Introduction to the LLL
 Algorithm and Its Applications" CRC Press, 2012.


# Examples
```jldoctest
julia> H= [1 2; 3 4];islllreduced(H)
false

julia> B,T=lll(H);islllreduced(B)
true

```
"""
function islllreduced(B)
# follows p 56 of Bremner
    Td = eltype(B)
    Xs,M = gso(B)
    isSizeReduced = all((M-I)[:]*2 .< one(Td))
    n = size(B,2)
    xiSq = zeros(float(Td),n)
    for i=1:n
        xiSq[i] = Xs[:,i]'*Xs[:,i]
    end
    deltas = zeros(float(Td),n-1)
    for i=2:n
        deltas[i-1] = xiSq[i]/xiSq[i-1]+M[i-1,i]^2
    end
    if isSizeReduced && minimum(deltas)>1/4
        return true
    else
        return false
    end
end

"""
    issizereduced(B)

Determine if the matrix B is size reduced or not.

# Examples
```jldoctest
julia> H= [1 2; 3 4];issizereduced(H)
false

julia> B,T = lll(H);issizereduced(B)
true

```
"""
function issizereduced(B)
    # to-do: replace gso with QR as on p 148 of LLL survey book
    Xs,M = gso(B)
    return all((M-I)[:]*2 .< one(eltype(B)))
end


"""
    Q,R = gso(H)

Find Gram-Schmidt orthogonal basis; H=Q*R
This is not quite a QR decomposition.

# Examples
```jldoctest
julia> H = [1 2; 3 4]; Q,R = gso(H)

```
"""
function gso(H::Matrix{Td}) where {Td}
# This is classic GS; see comparison between classic and modified in
# http://www.cis.upenn.edu/~cis610/Gram-Schmidt-Bjorck.pdf

# See also Fig 3.1 of "Lattice Basis Reduction: An Introduction to the LLL
# Algorithm and Its Applications" by Murray R. Bremner, CRC Press, 2012.

    Q = float.(H)            # Xs
    n = size(H,2)
    R = Matrix{float(Td)}(I,n,n) # M or μ

    for i=1:n
        for j=1:i-1
            R[j,i] = H[:,i]'*Q[:,j]/(Q[:,j]'*Q[:,j])
            Q[:,i]-= R[j,i]*Q[:,j]
        end
    end
    return Q,R
end

"""
    volume(B)

Volume of fundamental parallelepiped of a lattice with basis B.

# Examples
```jldoctest
julia> B = [1 2; 3 4]; volume(B)
1.9999999999999964

```
"""
volume(B) = sqrt(det(B'*B))
