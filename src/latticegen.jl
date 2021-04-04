
"""
    gen_qary!(matrix::Matrix{T},k::Int,q::Int) where {T}

...Maybe this should not try to operate in place...
"""
function gen_qary!(matrix::Matrix{T},k::Int,q) where {T}
    d,n = size(matrix)
    d==n || error("n!=d")

    matrix[1:d-k,1:d-k] = Matrix{T}(I,d-k,d-k)
    matrix[d-k+1:d,1:d-k] = rand(T(0):q-1,k,d-k)
    matrix[1:d-k,d-k+1:d] = zeros(T,k,d-k)
    matrix[d-k+1:d,d-k+1:d] = q*Matrix{T}(I,k,k)
end

"""
    b= gen_qary_b(T, d::Int,k::Int,b::Int)

Generate a q-ary lattice given an element type `T`, dimension `d`, parameter
`k`, and bit-depth `b`. Specifically, find a `d` by `d` matrix which has the
block structure `[I zeros(T,k,d-k); H q*I]]`, where the `k` by `d-k` matrix
H is sampled from `0:q-1` and q is sampled uniformly from `1:big(2)^b-1`

These bases correspond to the SIS/LWE q-ary lattices; see D. Micciancio and
O. Regev. Post-Quantum Cryptography. Chapter of Lattice-based Cryptography,
147-191 (2009) and latticegen in https://github.com/fplll/fplll

# Examples
```
julia> b=gen_qary_b(Int64,2,1,6)
2×2 Matrix{Int64}:
 1   0
 7  32

```
"""
function gen_qary_b(T, d::Int,k::Int,b::Int)
    if T!=BigInt && big(2)^b-1>typemax(T)
        error("Type $(T) can only handle b=$(log2(typemax(T)+1)) "*
              "bits, and you're asking for $b. Try again with "*
              "a different type or smaller bit depth b.")
    end
    q = T(rand(1:big(2)^b-1))
    matrix = Matrix{T}(undef,d,d)
    gen_qary!(matrix,k,q)
    return matrix
end
