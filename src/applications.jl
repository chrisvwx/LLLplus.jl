"""
    integerfeasibility(A,d,nullVecs=false)

Given a linear system `Ax=d`, return an integer vector `x` which satisfies the
system.  This is the *integer programming feasibility problem*.

If `nullVecs==true`, then as well as returning a solution `x`, also return a
matrix `xNull` of vectors in the null space of `A` which could be added to the
`x` vector to find a solution which satisfies a constraint such as `0 .≤ x
.≤ u`; see the paper below.

*Solving A System Of Diophantine Equations With Bounds On The Variables* by
Karen Aardal, Cor Hurkens, and Arjen Lenstra in Integer Programming and
Combinatorial Optimization, 6th International IPCO Conference, vol 1412, pp
229-242, 1998. See http://softlib.rice.edu/pub/CRPC-TRs/reports/CRPC-TR98782.pdf

# Examples
```jldoctest
julia> A=[10 1 -9; 1 8 8]; xtrue=[0; 2; 9]; d=A*xtrue;
julia> integerfeasibility(A,d)
3-element Array{Int64,1}:
 0
 2
 9

julia> n=20;m=30; A = rand(-10:10,n,m); xtrue = rand(0:10,m); d=A*xtrue;
julia> sum(abs.(xtrue-integerfeasibility(A,d)))
0

```
"""
function integerfeasibility(A::AbstractArray{Td,2},
                            d::AbstractArray{Td,1}, flag=false) where {Td<:Number}
    m,n=size(A)
    N1=10^3
    N2=10^4
    B = [[I           zeros(Td,n,1)];
         zeros(Td,1,n) N1;
         N2*A          -N2*d]

    Bhat,_ = lll(B)
    Bprime = Bhat[1:n+1,1:n-m+1]
    x_d = Bprime[1:n,end]

    if flag
        xNull=Bprime[1:n,1:n-m]
        return x_d,xNull
    else
        return x_d
    end
end



"""
    x = subsetsum(a,s)

For a vector of integers `a`, and an integer `s`, try to find a binary
vector `x` such that `x'*a=s`. We use the LLL algorithm to find the
solution.

This follows the technique described by Lagarias and Odlyzko  in 
"Solving Low-Density Subset Sum Problems"  in Journal of ACM, Jan 1985.
Code based on http://web.eecs.umich.edu/~cpeikert/lic15/lec05.pdf
The original technique came from Lagarias and Odlyzko.

It's odd that permuting the `a` vector in the second example given below
causes the alg to often not find a solution. The example doesn't have the
required density, so maybe we should be asking why it finds an answer at
all.

# Examples
```jldoctest
julia> a=[1,2,4,9,20,38]; s=30; x=subsetsum(a,s); s-x'*a
0

julia> a=[32771,65543,131101,262187,524387,1048759, # from Bremner p 117
          2097523,4195057,8390143,16780259,33560539,
          67121039,134242091,268484171,536968403];
julia> s=891221976; x=subsetsum(a,s); s-x'*a
0

julia> N=40;a=rand(1:2^BigInt(256),N);xtrue=rand(Bool,N); s=a'*xtrue; 
julia> setprecision(BigFloat,300); x=subsetsum(a,s); s-x'*a
0
```
"""
function subsetsum(a::AbstractArray{Ti,1},s::Ti) where {Ti<:Integer}
    # page numbers below refer to lecture note above

    n = length(a)
    
    flag = false;
    if s<sum(a)/2
        flag=true
        s = sum(a)-s
    end
    # b is the "B" parameter defined at bottom of page 2
    nt = Ti(n)
    b = ceil(sqrt(nt*2^nt))
    # BB is the 'B' matrix defined also at bottom of page 2
    BB = [Matrix{Ti}(I,n,n+1); -b*a' b*s]
    B,T = lll(BB);
    
    # CBP: just winging it from here...
    # Pick off the interesting parts of the B matrix
    Bp = B[1:n,1:n+1]
    # Check if any of the columns satisfy the requirement
    sums = abs.(a'*Bp)
    ixMatch = findall(sums.==s)
    # d = diag(Bp'*Bp)  # solution generally matches a small element of d

    binarySolution = false
    global xb
    for ix=1:length(ixMatch)
        # Pick off a column of B
        x = Bp[:,ixMatch[ix][2]]
        if !(any(x.<0)|any(x.>1))
            binarySolution= true
            xb = flag ? mod.(x .+ 1,2) : x
        end
        x = -x
        if !(any(x.<0)|any(x.>1))
            binarySolution= true
            xb = flag ? mod.(x .+ 1,2) : x
        end
    end

    if binarySolution
        return xb
    else
        if maximum(a)<2^(n^2/2)
            density = n/maximum(log2.(a))
            @printf("The density (%4.2f) of the 'a' vector is not as low as required\n",
                    density)
            @printf("(%4.2f) for Lagarias-Odlyzko to work. We'll look for a ",2/n)
            @printf("solution\nanyway... ")
        end

        if length(ixMatch)>=1
            print("A non-binary solution was found\n")
            return Bp[:,ixMatch[1][2]]
        else
            print("A solution was not found\n")
        end
    end
    return NaN
end
