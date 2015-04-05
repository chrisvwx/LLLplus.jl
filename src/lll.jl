function lll(H::Array{Complex{Float64},2},delta::FloatingPoint=3/4)
# (B,T,Q,R) = LLL(H,delta)
# Do Lenstra–Lenstra–Lovász lattice basis reduction of matrix H using
# parameter delta.  
#
# Input:     H     : input basis
#            delta : LLL reduction parameter
#
# Output:    B : LLL-reduced basis, B = H*T with QR decomposition B=Q*R
#            T : unimodular transformation matrix
#            Q : matrix with orthogonal columns of length one
#            R : upper triangular matrix
#
# Follows LLL in D. Wueben, et al, MMSE-Based Lattice-Reduction for Near-ML
# Detection of MIMO Systems International IEEE Workshop on Smart Antennas,
# Munich, March 2004.

if delta < 0.25 || delta > 1.0
    error("delta must be between 1/4 and 1.");
end

B = copy(H);
L = size(B,2);
(Q,R) = qr(B);
T = complex(eye(L));

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
        end
    end

    len = norm(R[lx-1:lx,lx]);
    if delta*abs(R[lx-1,lx-1])^2 > len^2
        
        # swap columns lx-1 and lx in B, T and R
        B[:,[lx-1,lx]]    = B[:,[lx,lx-1]];
        T[:,[lx-1,lx]]    = T[:,[lx,lx-1]];
        R[1:lx,[lx-1,lx]] = R[1:lx,[lx,lx-1]];

        # reconstruction of upper triangular structure by Givens rotation 
        # mutliplication with matrix Theta achieves R[lx,lx-1] = 0
        c     = R[lx-1,lx-1] / len;        # len = ||R[lx-1:lx,lx-1]|| after swapping
        s     = R[lx,lx-1]   / len;
        Theta = [c' s; -s c];
        
        R[lx-1:lx,lx-1:end] = Theta * R[lx-1:lx,lx-1:end];
        Q[:,lx-1:lx]       = Q[:,lx-1:lx] * Theta' ;
        lx                = max(lx-1,2);
    else
        lx = lx+1;
    end
end

return (B,T,Q,R)
end
