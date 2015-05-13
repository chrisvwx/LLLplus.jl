function hard_sphere{Td}(Y::Array{Td,2},H::Array{Td,2},Nc)
    # x=hard_sphere(y,H,Nc)
    #  Solve the problem argmin_x ||y-Hx||, where the elements of x are
    #  integers.  The input vector y is of dimension N, with H of dimension N
    #  by M, and the returned vector x of dimension M.
    # X=hard_sphere(Y,H)
    #  The input may be a matrix Y of dimension N by Ns, in which case
    #  the problem is solved for each column of Y, with the solutions in
    #  the columns of X. 
    #
    # Examples:
    #  X = hard_sphere([1 2]', [1 2; 3 4],2)
    (N,M) = size(H);

    if length(Nc)==1;
        Nc = ones(M,1)*Nc;
    end
    (Q,R) = qr(H);
    (N,Ns) = size(Y);
    Xh = zeros(Int,M,Ns);
    for ns = 1:Ns
        yp = (Q'*Y[:,ns] + R*ones(M,1).*(Nc-1))/2;
        xp = algII_smart(yp,R,Nc);
        #    Xh[:,ns] = xp*2 -(Nc-1); 
        Xh[:,ns] = xp
    end

    return Xh
end

#########################################################################%
function algII_smart(yp,R,Qc)
    # xh = algII_smart(yp,R,Qc)
    #
    # Find the closest xh for yp=R*xh + w, where w is some unknown noise and
    # x(i) is from the range [0,...,Qc(i)-1]. This implements the "Alg II-smart"
    # hard-decision sphere decoder as in Damen, El Gamal, and Caire, Trans It
    # 03.
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
                    ksi[ix-1] = (R[ix-1,ix:M]*xp[ix:M])[1];
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
