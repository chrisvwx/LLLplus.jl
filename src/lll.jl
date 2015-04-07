function lll{Td}(H::Array{Td,2},delta=3/4)
# (B,T,Q,R) = LLL(H,delta)
#  Do Lenstra–Lenstra–Lovász lattice basis reduction of matrix H using
#  parameter delta, which is optional with default delta=3/4.  The output is
#  B, an LLL-reduced basis; T, a unimodular transformation matrix such that
#  B= H*T; Q and R such that B=Q*R and Q is orthonormal, and R is upper
#  triangular.
#
# Example with large real matrix:
#   N=500;H = randn(N,N); include("lll.jl");(B,T,Q,R) = lll(H);
# Example with 2x2 complex matrix:
#   N=2;H = randn(N,N)+im*randn(N,N); include("lll.jl");(B,T,Q,R) = lll(H)

# By Christian Peel  (chris.peel@ieee.org)
# Last Modified: Sat 4 Apr 15, 9:26pm by peel

# Follows LLL in D. Wueben, et al, MMSE-Based Lattice-Reduction for Near-ML
# Detection of MIMO Systems International IEEE Workshop on Smart Antennas,
# Munich, March 2004.

if delta < 0.25 || delta > 1.0
    error("delta must be between 1/4 and 1.");
end

B = copy(H);
L = size(B,2);
(Q,R) = qr(B);
T = Td<:Complex? eye(Complex{Int},L): eye(Int,L);
    
lx  = 2;
while lx <= L
  
    # reduce lx-th column of B
    for k=lx-1:-1:1
        rk = R[k,lx]/R[k,k];
        mu = round(real(rk))+ im*round(imag(rk));
        if abs(mu)>0
            B[:,lx]   = B[:,lx]   - mu * B[:,k];
            R[1:k,lx] = R[1:k,lx] - mu * R[1:k,k];
            T[:,lx]   = T[:,lx]   - mu * T[:,k];
            # println("T=$(T)")
            # println("mu=$(mu)")
        end
    end

    nrm = norm(R[lx-1:lx,lx]);
    if delta*abs(R[lx-1,lx-1])^2 > nrm^2
        
        # swap columns lx-1 and lx in B, T and R
        B[:,[lx-1,lx]]    = B[:,[lx,lx-1]];
        T[:,[lx-1,lx]]    = T[:,[lx,lx-1]];
        R[1:lx,[lx-1,lx]] = R[1:lx,[lx,lx-1]];

        # upper triangular by Givens rotation 
        # mutl with matrix Theta achieves R[lx,lx-1] = 0
        cc = R[lx-1,lx-1] / nrm; # nrm = ||R[lx-1:lx,lx-1]|| after swapping
        ss = R[lx,lx-1]   / nrm;
        Theta = [cc' ss; -ss cc];
        
        R[lx-1:lx,lx-1:end] = Theta * R[lx-1:lx,lx-1:end];
        Q[:,lx-1:lx]       = Q[:,lx-1:lx] * Theta' ;
        lx                = max(lx-1,2);
    else
        lx = lx+1;
    end
end

return (B,T,Q,R)
end
