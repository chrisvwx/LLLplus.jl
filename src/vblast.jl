function vblast{T}(H::Array{T,2},mu=Inf)
# (W,P,B) = vblast(H)
#  Returns a VBLAST matrix decomposition of H such that 
#  H = pinv(W)*B*P' or B = W*H*P.  Here P is a permutation matrix, B is
#  lower triangular with ones on the diagonal, and W has orthogonal rows.
# (W,P,B) = vblast(H,mu)
#  If an SNR argument mu is passed in, a regularized ("MMSE")
#  decomposition is done, with the result that W will no longer have
#  orthogonal rows and B is no longer lower triangular.
#
# Example: For decoding with N receivers
#  N=4;H = rand(N,N); include("vblast.jl");(W,P,B) = vblast(H)
# Example: For MMSE decoding with N receivers mu=rho/N 
#  (W,P,B) = vblast(H,rho/size(H,1)); 

(K,N) = size(H);

W = zeros(T,K,K);
P = zeros(Int,N,N);
Hp = H;
k = zeros(Int,K,1);
indexH = [1:N;];
for i=1:N
    G = inv(Hp'*Hp + eye(N-i+1)/mu)*Hp'; # Do pseudoinverse
    gn = diag(G*G');
    (mx,mi) = findmin(abs(gn))
    k[i] = indexH[mi];
    W[i,:] = G[mi,:];
    P[k[i],i]=1;
    indexH = setdiff(indexH,k[i]);
    Hp = H[:,indexH];
end
B = W*H*P  # Calculating B may not be required in some real-world scenarios
return (W,P,B)
end
