function seysen{Td}(H::Array{Td,2})
# [B,T,B_dual,num_it] = seysen{Td}(H::Array{Td,2})
#
# Do greedy Seysen lattice reduction on the matrix H (real or
# complex-valued), returning T a unimodular matrix that reduces H;
# B, the reduced lattice basis (i.e. B = H*T);
# B_dual,  dual lattice basis (i.e., B_dual = pinv(B)); and
# the number of iterations (the number of basis updates).

# Follows LLL in D. Wueben, et al, MMSE-Based Lattice-Reduction for Near-ML
# Detection of MIMO Systems International IEEE Workshop on Smart Antennas,
# Munich, March 2004.

if Td<:BigInt || Td<:BigFloat
    Ti= Td<:Complex?Complex{BigInt}:BigInt
else
    Ti= Td<:Complex?Complex{Int}:Int
end
roundf(r) = Td<:Complex?round(real(r)) + im*round(imag(r)):round(r);

#H = H*1.0

# get size of input
(n, m) = size(H); # m-dimensional lattice in an n-dimensional space

# initialization, outputs
B      = copy(H);    # reduced lattice basis
num_it = 0;          # number of iterations
T      = eye(Ti,m);  # unimodular matrix
A      = (H'*H)*1.0;       # Gram matrix of H
Adual  = inv(A);     # Inverse gram matrix of H
B_dual = H*Adual;    # Dual basis

# calculate all possible update values Λ[s,t)
# and their corresponding reduction Δ[s,t) in Seysen's measure
Λ = zeros(Ti,m,m);
Δ = zeros(Td,m,m)*0.0;

for s = 1:m
    for t = 1:m
        if s != t
	    x = 0.5*(Adual[t,s]/Adual[s,s]-A[t,s]/A[t,t]);
            Λ[s,t] = roundf(x);
	    AbsΛ   = abs(Λ[s,t])^2;
	    if AbsΛ != 0
	        zw = real(Λ[s,t])*real(x)+imag(Λ[s,t])*imag(x);
	        Δ[s,t]  = Adual[s,s]*A[t,t]*(2*zw-AbsΛ);
	    end
        end
    end
end # - end calculation of Λ and Δ

# find maximum reduction in Seysen's measure (greedy approach)
(zw,max_ind) = findmax(abs(Δ[:]));
(s, t)        = ind2sub((m,m),max_ind);

# init loop
do_reduction = true;

# if no improvement can be achieved
if Δ[s,t] == 0
    do_reduction = false;
end

# Lattice reduction loop of subsequent basis updates
while do_reduction

    num_it = num_it + 1;
    
    # perform basis update
    B[:,s] = B[:,s] + Λ[s,t]*B[:,t];

    # compute corresponding unimodular trasformation matrix
    T[:,s] = T[:,s] + Λ[s,t]*T[:,t];  # basis transformation matrix

    # update corresponding dual basis
    B_dual[:,t] = B_dual[:,t] - Λ[s,t]'*B_dual[:,s];

    # Update Gram and inverse Gram matrix
    for ind = 1:m

        # update Gram matrix
        if ind != s
            A[s,ind] = A[s,ind]+Λ[s,t]'*A[t,ind];
            A[ind,s] = A[s,ind]';
        else
            A[s,ind] = norm(B[:,s]).^2;
        end

        # update inverse Gram matrix
        if ind != t
            Adual[t,ind] = Adual[t,ind]-Λ[s,t]*Adual[s,ind];
            Adual[ind,t] = Adual[t,ind]';
        else
            Adual[t,ind] = norm(B_dual[:,t])^2;
        end

    end # - end update Gram und inverse Gram matrix
    
    # update all possible update values Λ[s,t]
    # and their corresponding reduction Δ[s,t] in Seysen's measure
    for ind1 = 1:m
        for ind2 = 1:m
            if ((ind1==s) | (ind1==t) | (ind2==s) | (ind2==t)) & (ind1!=ind2)
                x = 0.5*(Adual[ind2,ind1]/Adual[ind1,ind1]-A[ind2,ind1]/
                         A[ind2,ind2]);
                try
                    Λ[ind1,ind2] = roundf(x);
                catch
                    println("x = $(x), typeof(x) = $(typeof(x))")
                    println("A[ind2,ind2] = $(A[ind2,ind2]), typeof(A[ind2,ind2]) = $(typeof(A[ind2,ind2]))")
                    println("roundf(x)) = $(roundf(x))"*
                            "typeof(roundf(x))) = $(typeof(roundf(x)))")
                    println("Λ[ind1,ind2]) = $(Λ[ind1,ind2])"*
                            "typeof(Λ[ind1,ind2])) = $(typeof(Λ[ind1,ind2]))")
                    println(" ")
                end
                AbsΛ = abs(Λ[ind1,ind2])^2;
                if AbsΛ != 0
                    zw = real(Λ[ind1,ind2])*real(x)+imag(Λ[ind1,ind2])*imag(x);
                    Δ[ind1,ind2] = Adual[ind1,ind1]*A[ind2,ind2]*(2*zw-AbsΛ);
                else
                    Δ[ind1,ind2] = 0;
                end
            end
        end
    end # - end update Λ and Δ

    # find maximum reduction in Seysen's measure (greedy approach)
    (zw, max_ind) = findmax(abs(Δ[:]));
    (s, t)        = ind2sub((m,m),max_ind);

    # if no reduction is possible, exit loop
    if Δ[s,t] == 0
        do_reduction = false;
    end
	
end # - end lattice reduction loop


return (B,T,B_dual,num_it)
end
