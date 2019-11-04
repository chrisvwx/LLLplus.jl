"""
    B,T,B_dual,num_it = seysen(H::Array{Td,2}) where Td

Do greedy  Seysen lattice reduction  on the  matrix `H`, returning  `B`, the
reduced lattice basis;  `T` a unimodular matrix that reduces  `H` (i.e. `B =
H*T`); `B_dual`, dual lattice basis (i.e., `B_dual = pinv(B)`); and num_it the number
of iterations (basis updates). See also [`lll`](@ref).

Follows Seysen algorithm in "Lattice Reduction - A Survey with Applications
in Wireless Communications" by D. Wuebben, et al, IEEE Signal Processing
Magazine, 2011.

# Examples
```jldoctest
julia> H= [1 2; 3 4];B,T = seysen(H); B
2×2 Array{Int64,2}:
 -1  1
  1  1

julia> H= BigFloat.([1.5 2; 3 4]) .+ 2im; B,_= seysen(H); B
2×2 Array{Complex{BigFloat},2}:
 0.0+1.0im  0.50+0.0im
 0.0+0.0im   1.0+0.0im

```
"""
function seysen(H::Array{Td,2}) where {Td}

if Td<:BigInt || Td<:BigFloat
    Ti= Td<:Complex ? Complex{BigInt} : BigInt
else
    if Td<:Integer
        Ti= Td<:Complex ? Complex{Td} : Td
    else
        Ti= Td<:Complex ? Complex{Int} : Int
    end
end

n,m = size(H); # m-dimensional lattice in an n-dimensional space

# initialization, outputs
B      = copy(H);    # reduced lattice basis
num_it = 0;          # number of iterations
T      = Matrix{Ti}(I, m, m)  # unimodular matrix
A      = (H'*H)*1.0;          # Gram matrix of H
Adual  = inv(A);     # Inverse gram matrix of H
B_dual = H*Adual;    # Dual basis

# calculate update values Λ[s,t] and the corresponding reduction Δ[s,t] in
# Seysen's measure 
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
zw,max_ind = findmax(abs.(Δ[:]));
s,t        = Tuple(CartesianIndices((m,m))[max_ind])

# init loop
do_reduction = true;

if Δ[s,t] == 0  # if no improvement can be achieved
    do_reduction = false;
end

# Lattice reduction loop of subsequent basis updates
while do_reduction

    num_it = num_it + 1;

    
    B[:,s]  = B[:,s] + Λ[s,t]*B[:,t];    # perform basis update
    T[:,s]  = T[:,s] + Λ[s,t]*T[:,t];    # updater unimodular transformation
    B_dual[:,t] = B_dual[:,t] - Λ[s,t]'*B_dual[:,s];  # update dual basis

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
    end

    # update all possible update values Λ[s,t]
    # and their corresponding reduction Δ[s,t] in Seysen's measure
    for ind1 = 1:m
        for ind2 = 1:m
            if ((ind1==s) | (ind1==t) | (ind2==s) | (ind2==t)) & (ind1!=ind2)
                x = 0.5*(Adual[ind2,ind1]/Adual[ind1,ind1]-A[ind2,ind1]/
                         A[ind2,ind2]);

                # If you get an "InexactError:..." here, try a float type
                # (or float+int type) with more precision. For some reason
                # Seysen is more sensitive than LLL to bit precision.
                Λ[ind1,ind2] = Ti(roundf(x));

                AbsΛ = abs(Λ[ind1,ind2])^2;
                if AbsΛ != 0
                    zw = real(Λ[ind1,ind2])*real(x)+imag(Λ[ind1,ind2])*imag(x);
                    Δ[ind1,ind2] = Adual[ind1,ind1]*A[ind2,ind2]*(2*zw-AbsΛ);
                else
                    Δ[ind1,ind2] = 0;
                end
            end
        end
    end

    # find maximum reduction in Seysen's measure (greedy approach)
    zw,max_ind = findmax(abs.(Δ[:]));
    s,t = Tuple(CartesianIndices((m,m))[max_ind])

    # if no reduction is possible, exit loop
    if Δ[s,t] == 0
        do_reduction = false;
    end

end # while do_reduction


return B,T,B_dual,num_it
end
