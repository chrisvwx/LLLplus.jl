"""
    x=hardsphere(y,H,Nc)

  Solve the problem argmin_x ||y-Hx||.  The input vector y is of
  length N, with H of dimension N by M, and the returned vector x of
  length M.  If Nc is even the elements of x are integers from
  [0:Nc-1]*2-(Nc-1); odd Nc is not yet supported. For Nc=2, we're
  searching the 2PAM constellation of [-1,1], and for nInt=4 we're
  searching [-3,-1,1,3].

  Digital communication researchers often call the Fincke-Pohst search
  algorithm used in this function as the "sphere-decoder" because of its use
  in decoding digital signals and for the use of a multi-dimensional sphere
  to constrain search complexity. This function uses only hard-decision
  decoding, hence the term "hardsphere".

X=hardsphere(Y,H,Nc)

  The input may be a matrix Y of dimension N by Ns, in which case the
  problem is solved for each column of Y, with the solutions in the
  columns of X.

# Examples:
```jldoctest
julia> X = hardsphere([1. 1]', [1. 2; 3 4],2)
2×1 Matrix{Int64}:
 -1
  1

```
"""
function hardsphere(Y::AbstractArray{Td,2},H::AbstractArray{Td,2},Nc::Integer) where {Td}
    (N,M) = size(H);

    Qc = ones(M)*Nc;
    if iseven(Nc)
        xoffset = Qc.-1
        xmult = 2
    else
        @error("Odd Nc is not yet supported.")
        # To find a failing case, remove @error above, and compare w full
        # search. The fix is likely minor (once it is found).
        xoffset = floor.(Qc/2)
        xmult = 1
    end
    Q,R = qr(H);
    N,Ns = size(Y);
    Xh = zeros(Int,M,Ns);
    for ns = 1:Ns
        yp = (Q'*Y[:,ns] + R*ones(M,1).*xoffset)/xmult;
        xp = algIIsmart(yp,R,Qc);
        Xh[:,ns] = xp*xmult - xoffset; 
    end

    return Xh
end

"""
    hard_sphere(...) = hardsphere(...)

See hardsphere; in a future version hard_sphere will be removed.
"""
hard_sphere(x...) = hardsphere(x...)

"""
    xh = algIIsmart(yp,R,Qc)

Find the minimizing xh = argmin_x ||yp- R*x||, where `yp` is a length N
vector, `R` is NxM, `xh` will be a length M vector, and x(i) is
from the range [0,...,Qc(i)-1]. R is assumed to be upper triangular.

This implements the "Alg II-smart" hard-decision sphere decoder as in "On
Maximum-Likelihood Detection and the Search for the Closest Lattice Point,"
M. O. Damen, H. El Gamal, and G. Caire, Transactions Information Theory,
v49, Oct 2003. This function is not exported from LLLplus at present and may
disappear in the future.

# Examples:
```jldoctest
julia> X = LLLplus.algIIsmart([5. 8]', [1. 2; 0 4],[4,4])
2×1 Matrix{Float64}:
 1.0
 2.0

```
"""
function algIIsmart(yp,R,Qc)
    (N,M) = size(R);

    xh = zeros(M,1); # the best solution
    xp = zeros(M,1); # work vector (a possible solution)
    ksi = zeros(M,1);
    TT = zeros(M,1);
    delta = zeros(M,1);

    # Step 1: initialization
    ix = M;
    TT[ix] = 0;
    ksi[ix] = 0;
    dc = Inf;
    state = 2;
    terminate= false;

    # show(Qc)
    # println(" ")
    #    @printf("state=%d, ix=%d, xp[ix] = %6.3f, xh[ix] = %6.3f, "*
    #       "delta[ix] = %6.3f\n",state,ix,xp[ix],xh[ix],delta[ix])
    while ~terminate

        if state==2 # Step 2: DFE
            xp[ix] = round((yp[ix]-ksi[ix])/R[ix,ix]);
            delta[ix] = sign(yp[ix]-ksi[ix] - R[ix,ix]*xp[ix]);
            state = 3;
        elseif state== 3 # Step 3: Main loop
            ytmp = yp[ix]-ksi[ix] - R[ix,ix]*xp[ix];
            # @printf("ytmp'*ytmp = %6.3f, TT[ix] = %6.3f, "*
            #   "TT[ix] + ytmp'*ytmp = %6.3f, dc=%6.3f\n",
            #   ytmp'*ytmp,TT[ix],TT[ix] + ytmp'*ytmp,dc)
            if TT[ix] + ytmp'*ytmp >= dc
                # We're outside the sphere
                state = 4;
                # elseif !(0<=real(xp[ix])<Qc[ix]) || !(0<=imag(xp[ix])<Qc[ix])
            elseif !(0 <= xp[ix] < Qc[ix])
                # inside sphere but outside signal set boundary
                state = 6;
            else # inside sphere and inside signal set
                if ix>1
                    # ksi[ix-1] = (conj(R[ix-1,ix:M])*xp[ix:M])[1];
                    ksi[ix-1] = (R[ix-1,ix:M]'*xp[ix:M])[1];
                    TT[ix-1] = TT[ix] + ytmp'*ytmp;
                    ix = ix-1;
                    state = 2;
                else
                    state = 5;
                end
            end
        elseif state==4 # Step 4: Check termination criteria
            if ix==M
                terminate=true;
            else
                ix = ix+1;
                state = 6;
            end
        elseif state==5 # Step 5: A valid point is found
            ytmp = yp[1] -ksi[1] - R[1,1]*xp[1];
            dc = TT[1] + ytmp'*ytmp;
            xh[:] = xp[:];
            ix = ix+1;
            state = 6;
        elseif state==6 # Step 6: Schnorr-Euchner
            xp[ix] = xp[ix] + delta[ix];
            delta[ix] = -delta[ix] -sign(delta[ix]);
            state = 3;
        else
            error("Unknown state");
        end
        #@printf("state=%d, ix=%d, xp = [%1.0f %1.0f %1.0f %1.0f], "*
        #    xh = [%1.0f %1.0f %1.0f %1.0f], δ[ix] = %2.0f\n",
        #    state,ix,xp[1],xp[2],xp[3],xp[4],xh[1],xh[2],xh[3],xh[4],delta[ix])
        #    @printf("state=%d, ix=%d, xp = [%1.0f %1.0f], xh = [%1.0f "*
        #      "%1.0f],  δ[ix] = %2.0f\n",state,ix,xp[1],xp[2],xh[1],xh[2],delta[ix])
    end
    #println("xh=$(xh)")
    return xh
end
