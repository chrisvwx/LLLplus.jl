function seysen{Td}(H::Array{Td,2})
# [T,H_red,H_red_dual,num_it] = seysen{Td}(H::Array{Td,2})
#
# Do greedy Seysen lattice reduction on the matrix H (real or
# complex-valued), returning T a unimodular matrix that reduces H;
# H_red, the reduced lattice basis (i.e. H_red = H*T);
# H_red_dual,  dual lattice basis (i.e., H_red_dual = pinv(H_red)); and
# the number of iterations (the number of basis updates).

# By Chris Peel  (chris.peel@ieee.org)
# Last Modified: Thu 16 Apr 15, 4:40pm

# Follows LLL in D. Wueben, et al, MMSE-Based Lattice-Reduction for Near-ML
# Detection of MIMO Systems International IEEE Workshop on Smart Antennas,
# Munich, March 2004.

# get size of input
(n, m) = size(H); # m-dimensional lattice in an n-dimensional space

if isempty(H)
    error("Input basis is empty");
elseif rank(H) < m
    error("Input basis has not full column rank");
end
Ti = Td<:Complex? Complex{Int}: Int; # Integer type: Complex or not
Tf = Td<:Complex? Complex{FloatingPoint}: FloatingPoint;
H = convert(Array{Tf,2},H);

# initialization, outputs
H_red      = copy(H);    # reduced lattice basis
num_it     = 0;          # number of iterations
T          = eye(Ti,m);  # unimodular matrix
A          = H'*H;       # Gram matrix of H
Adual      = inv(A);     # Inverse gram matrix of H
H_red_dual = H*Adual;    # Dual basis

# calculate all possible update values Λ[s,t)
# and their corresponding reduction Δ[s,t) in Seysen's measure
Λ  = zeros(Ti,m,m);
Δ   = zeros(Tf,m,m);

for s = 1:m
    for t = 1:m
        if s != t
	    x = 0.5*(Adual[t,s]/Adual[s,s]-A[t,s]/A[t,t]);
            Λ[s,t] = round(real(x))+ im*round(imag(x));
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
    H_red[:,s] = H_red[:,s] + Λ[s,t]*H_red[:,t];

    # compute corresponding unimodular trasformation matrix
    T[:,s] = T[:,s] + Λ[s,t]*T[:,t];  # basis transformation matrix

    # update corresponding dual basis
    H_red_dual[:,t] = H_red_dual[:,t] - Λ[s,t]'*H_red_dual[:,s];

    # Update Gram and inverse Gram matrix
    for ind = 1:m

        # update Gram matrix
        if ind != s
            A[s,ind] = A[s,ind]+Λ[s,t]'*A[t,ind];
            A[ind,s] = A[s,ind]';
        else
            A[s,ind] = norm(H_red[:,s]).^2;
        end

        # update inverse Gram matrix
        if ind != t
            Adual[t,ind] = Adual[t,ind]-Λ[s,t]*Adual[s,ind];
            Adual[ind,t] = Adual[t,ind]';
        else
            Adual[t,ind] = norm(H_red_dual[:,t])^2;
        end

    end # - end update Gram und inverse Gram matrix
    
    # update all possible update values Λ[s,t]
    # and their corresponding reduction Δ[s,t] in Seysen's measure
    for ind1 = 1:m
        for ind2 = 1:m
            if ((ind1==s) | (ind1==t) | (ind2==s) | (ind2==t)) & (ind1!=ind2)
                x = 0.5*(Adual[ind2,ind1]/Adual[ind1,ind1]-A[ind2,ind1]/
                         A[ind2,ind2]);
                #            Λ[ind1,ind2] = round(x);
                Λ[ind1,ind2] = round(real(x))+ im*round(imag(x));
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


return (T,H_red,H_red_dual,num_it)
end
