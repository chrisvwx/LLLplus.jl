function seysen{Td}(H::Array{Td,2})
# [T,H_red,H_red_dual,num_it] = Seysen(H)
# This function performs (greedy) Seysen lattice reduction (see [1]) 
# 
# INPUT: 
# H             ... input basis (in general, complex-valued) 
#
# OUTPUT:  
#
# T             ... unimodular transformation matrix reducing H  
# H_red         ... reduced lattice basis, i.e., H_red = H*T 
# H_red_dual    ... dual lattice basis, i.e., H_red_dual = pinv(H_red)   
# num_it        ... number of iterations (number of basis updates)  

# See D. Wuebben, D. Seethaler, J. Jalden, and G. Matz, 
# "Lattice Reduction - A Survey with Applications in Wireless Communications" 
# IEEE Signal Processing Magazine, March 2011

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
H_red           = copy(H);    # reduced lattice basis 
num_it          = 0;          # number of iterations 
T               = eye(Ti,m);  # unimodular matrix  
A               = H'*H;       # Gram matrix of H 
Adual           = inv(A);     # Inverse gram matrix of H 
H_red_dual      = H*Adual;    # Dual basis 


# calculate all possible update values Lambda[s,t) 
# and their corresponding reduction Delta[s,t) in Seysen's measure 
Lambda  = zeros(Ti,m,m);
Delta   = zeros(Tf,m,m);
    
for s = 1:m 
    for t = 1:m 
        if s != t
	    x = 0.5*(Adual[t,s]/Adual[s,s]-A[t,s]/A[t,t]); 
# println("s = $(s), typeof(s) = $(typeof(s))")
# println("t = $(t), typeof(t) = $(typeof(t))")
# println("Lambda = $(Lambda), typeof(Lambda) = $(typeof(Lambda))")
            Lambda[s,t] = round(real(x))+ im*round(imag(x));
	    AbsLambda   = abs(Lambda[s,t])^2;
	    if AbsLambda != 0
	        zw = real(Lambda[s,t])*real(x)+imag(Lambda[s,t])*imag(x);
# println("typeof(Adual[s,s]) = $(typeof(Adual[s,s]))")
# println("typeof(A[t,t]) = $(typeof(A[t,t]))")
# println("zw = $(zw),  AbsLambda = $(AbsLambda)")
# println("Delta[s,t] = $(Delta[s,t]), typeof(Delta[s,t]) = $(typeof(Delta[s,t]))")
# println("expr = $(Adual[s,s]*A[t,t]*(2*zw-AbsLambda))")
# println("expr = $(typeof(Adual[s,s]*A[t,t]*(2*zw-AbsLambda)))")
	        Delta[s,t]  = Adual[s,s]*A[t,t]*(2*zw-AbsLambda);
	    end
        end
    end
end # - end calculation of Lambda and Delta

# find maximum reduction in Seysen's measure (greedy approach)
# println("Delta[:] = $(Delta[:])")
# println("abs(Delta[:]) = $(abs(Delta[:]))")
(zw,max_ind) = findmax(abs(Delta[:]));
(s, t)        = ind2sub((m,m),max_ind);

# init loop
do_reduction = true;

# if no improvement can be achieved
if Delta[s,t] == 0 
   do_reduction = false; 
end

# Lattice reduction loop of subsequent basis updates
while do_reduction 

   # increase iteration counter
   num_it = num_it + 1;
	
   # perform basis update
   H_red[:,s] = H_red[:,s] + Lambda[s,t]*H_red[:,t];			

   # compute corresponding unimodular trasformation matrix
   T[:,s] = T[:,s] + Lambda[s,t]*T[:,t];  # basis transformation matrix 

   # update corresponding dual basis
   H_red_dual[:,t] = H_red_dual[:,t] - Lambda[s,t]'*H_red_dual[:,s];
   
   # Update Gram and inverse Gram matrix
   for ind = 1:m 
   
      # update Gram matrix
      if ind != s
         A[s,ind] = A[s,ind]+Lambda[s,t]'*A[t,ind]; 	   
         A[ind,s] = A[s,ind]';
      else 
         A[s,ind] = norm(H_red[:,s]).^2;
      end
   
      # update inverse Gram matrix
      if ind != t
         Adual[t,ind] = Adual[t,ind]-Lambda[s,t]*Adual[s,ind]; 
         Adual[ind,t] = Adual[t,ind]';	   
      else
         Adual[t,ind] = norm(H_red_dual[:,t])^2;
      end
	  
   end # - end update Gram und inverse Gram matrix
						
   # update all possible update values Lambda[s,t] 
   # and their corresponding reduction Delta[s,t] in Seysen's measure 
   for ind1 = 1:m 
      for ind2 = 1:m 
         if ((ind1 == s) | (ind1 == t) | (ind2 == s) | (ind2 == t)) & (ind1 != ind2) 
            x = 0.5*(Adual[ind2,ind1]/Adual[ind1,ind1]-A[ind2,ind1]/A[ind2,ind2]); 
#            Lambda[ind1,ind2] = round(x);
            Lambda[ind1,ind2] = round(real(x))+ im*round(imag(x));
            AbsLambda = abs(Lambda[ind1,ind2])^2;
            if AbsLambda != 0
               zw = real(Lambda[ind1,ind2])*real(x)+imag(Lambda[ind1,ind2])*imag(x);
               Delta[ind1,ind2] = Adual[ind1,ind1]*A[ind2,ind2]*(2*zw-AbsLambda);
            else
               Delta[ind1,ind2] = 0;
            end
         end
      end
   end # - end update Lambda and Delta
    
   # find maximum reduction in Seysen's measure (greedy approach)
   (zw, max_ind) = findmax(abs(Delta[:]));
   (s, t)        = ind2sub((m,m),max_ind); 
    
   # if no reduction is possible, exit loop
   if Delta[s,t] == 0 
      do_reduction = false;
   end
	
end # - end lattice reduction loop


return (T,H_red,H_red_dual,num_it)
end
