"""
    W,P,B = vblast(H)

Find a VBLAST decomposition of `H` such that `H = pinv(W)*B*P'` or `B =
W*H*P`.  Here `P` is a permutation matrix, `B` is lower triangular with ones
on the diagonal, and `W` has orthogonal rows.

    W,P,B = vblast(H,mu)

If an SNR argument `mu` is passed in, a regularized ("MMSE") decomposition
is done, with the result that `W` will no longer have orthogonal rows and
`B` is no longer lower triangular.


# Examples
```jldoctest
julia> H= [1. 2; 3 4];W,_ = vblast(H); W
2×2 Matrix{Float64}:
 1.5  -0.5
 0.1   0.3

julia> H= BigFloat.([1.5 2; 3 4]) .+ 2im; W,_= vblast(H); W
2×2 Matrix{Complex{BigFloat}}:
      -2.0+3.0im            2.0-1.5im     
 0.0779221-0.103896im  0.155844-0.103896im

```
"""
function vblast(H::Array{T,2},mu=Inf) where {T}

(K,N) = size(H);

W = zeros(T,K,K);
P = zeros(Int,N,N);
Hp = H;
k = zeros(Int,K,1);
indexH = [1:N;];
for i=1:N
    G = inv(Hp'*Hp + I/mu)*Hp'; # Do pseudoinverse
    gn = diag(G*G');
    (mx,mi) = findmin(abs.(gn))
    k[i] = indexH[mi];
    W[i,:] = G[mi,:];
    P[k[i],i]=1;
    indexH = setdiff(indexH,k[i]);
    Hp = H[:,indexH];
end
B = W*H*P  # Calculating B may not be required in some real-world scenarios
return (W,P,B)
end
