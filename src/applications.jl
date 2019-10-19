"""
    integerfeasibility(A,d,nullVecs=false)

Given a linear system `Ax=d`, return an integer vector `x` which satisfies the
system.  This is the *integer programming feasibility problem*.

If `nullVecs==true`, then as well as returning a solution `x`, also return a
matrix `xNull` of vectors in the null space of `A` which could be added to the
`x` vector to find a solution which satisfies a constraint such as `0 .≤ x
.≤ u`; see the paper below.  

This is not a robust tool, just a demo.

"Solving A System Of Diophantine Equations With Bounds On The Variables" by
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

julia> sum(abs.(xtrue - integerfeasibility(A,d) ))
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
solution. This is not a robust tool, just a demo.

This follows the technique described by Lagarias and Odlyzko  in 
"Solving Low-Density Subset Sum Problems"  in Journal of ACM, Jan 1985.
Code based on http://web.eecs.umich.edu/~cpeikert/lic15/lec05.pdf
We can likely get better results using techniques described and referenced in
https://www-almasty.lip6.fr/~joux/pages/papers/ToolBox.pdf

It's odd that permuting the `a` vector in the second example given below
causes the alg to often not find a solution. The example doesn't have the
required density, so maybe we should be asking why it finds an answer at
all.

# Examples
```jldoctest
julia> a=[1,2,4,9,20,38]; s=30; x=subsetsum(a,s); s-x'*a
0.0

julia> a=[32771,65543,131101,262187,524387,1048759, # from Bremner p 117
          2097523,4195057,8390143,16780259,33560539,
          67121039,134242091,268484171,536968403];

julia> s=891221976; x=subsetsum(a,s); s-x'*a
0.0

julia> N=40;a=rand(1:2^BigInt(256),N);xtrue=rand(Bool,N); s=a'*xtrue; 

julia> setprecision(BigFloat,300); x=subsetsum(a,s); s-x'*a
0.0

```
"""
function subsetsum(a::AbstractArray{Ti,1},ss::Ti,returnBinary=false) where {Ti<:Integer}
    # page numbers below refer to lecture note above

    n = length(a)
    
    if ss<sum(a)/2
        flag=true
        s = sum(a)-ss
    else
        flag = false
        s = ss
    end
    # b is the "B" parameter defined at bottom of page 2
    nt = Ti(n)
    b = ceil(sqrt(nt*2^nt))
    # BB is the 'B' matrix defined also at bottom of page 2
    BB = [Matrix{Ti}(I,n,n+1); -b*a' b*s]
    B,T = lll(BB)

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
        x = mdsubsetsum(a,ss,.5,1)
        if ~ismissing(x)
            print("A solution was found via mdsubsetsum\n")
            return x
        end

        if maximum(a)<2^(n^2/2)
            density = n/maximum(log2.(a))
            @printf("The density (%4.2f) of the 'a' vector is not as low as required\n",
                    density)
            @printf("(%4.2f) for Lagarias-Odlyzko to work. \n",2/n)
        end

        if returnBinary && length(ixMatch)>=1
            print("A non-binary solution was found; check that it's correct\n")
            return Bp[:,ixMatch[1][2]]
        else
            print("A solution was not found\n")
        end
    end
    return missing.*Bp[:,1]
end


"""
    x = mdsubsetsum(a,sM,ratio=.5,Kpm=3)

For a vector of integers `a`, and an integer `sM`, try to find a binary
vector `x` such that `x'*a=s` using the technique from "Multidimensional
subset sum problem" [1][2]. A major goal of the technique is to solve
problems in there are about 50% ones in `x`; other ratios of ones to zeros
can be specified in `ratio`.  The thesis also suggests searching `Kpm=3`
values around the nominal k. This technique is related to that in
[`subsetsum`](@ref) in that both use the LLL algorithm.  This is not a
robust tool, just a demo.

[1] https://scholarworks.rit.edu/theses/64/
[2] https://pdfs.semanticscholar.org/21a7/c2f9ff29507f1153aefcca04d1cd308e45c0.pdf

# Examples
```jldoctest
julia> a=[1,2,4,9,20,38]; s=30; x=mdsubsetsum(a,s); s-x'*a
0.0

julia> a=[32771,65543,131101,262187,524387,1048759, # from Bremner p 117
          2097523,4195057,8390143,16780259,33560539,
          67121039,134242091,268484171,536968403];

julia> sM=891221976; x=mdsubsetsum(a,sM); sM-x'*a

julia> N=40;a=rand(1:2^BigInt(256),N);xtrue=rand(Bool,N); s=a'*xtrue;

julia> setprecision(BigFloat,300); x=mdsubsetsum(a,s); s-x'*a
```
"""
function mdsubsetsum(a::AbstractArray{Ti,1},sM::Ti,ratio=.5,Kpm=3) where {Ti<:Integer}

    n = Ti(length(a))

    # r is heuristic value for r from end of Section 11 of thesis above
    r = Ti(2^round(.65*log2(maximum(a))))
    c = 2^10 # a heuristic; look right below eq (10)

    p = a .÷ r # See start of Section 10
    s = a-p.*r
    k0 = Ti(round(sum(s)*ratio/r)) # Equation (9) with ratio=.5 and j=1.
    p0 = sM ÷ r
    m = sM - p0*r

    for k = k0-Kpm:k0+Kpm
        # Equation (14)
        Bt = [Matrix{Ti}(I,n,n)*2 c*s       c*p      zeros(Ti,n);
              ones(Ti,1,n)        c*(k*r+m) c*(p0-k) 1]
        BB = Bt'  # Notation in paper is row-based; my brain is column-based
        B,T = lll(BB)

        for nx = 1:n
            if abs(B[n+3,nx])==1 && B[n+2,nx]==0 && B[n+1,nx]==0 &&
                      all((B[1:n,nx] .==1) .| (B[1:n,nx] .==-1))
                xhat = abs.(B[1:n,nx] .- B[n+3,nx])/2
                return Ti.(xhat)
            end
        end
    end
    @warn "Solution not found"
    return missing
end


"""
    rationalapprox(x::AbstractArray{<:Real,1},M,Ti=BigInt,verbose=false)

For a vector of Reals `x`, and an integer `M`, find an integer q such that
`maximum(abs.(x*q-round.(x*q)))` is small; the vector `x` is approximated by
`round.(x*q)//q`.  The integer `q` is less than or equal to `M` and the
approximation satisfies `max(abs.(x*q-round.(x*q)))≤sqrt(5)*2^(n/4 -
5)*M^(-1/n)`; this equation comes from the paper below.  The LLL algorithm
reduction is used to find the solution. The approximation vector is returned.
This is also known as "simultaneous diophantine approximation"; see for
example the title of the Hanrot paper below.

This is not a robust tool, just a demo.

"LLL: A Tool for Effective Diophantine Approximation" by Guillaume Hanrot in
the book "The LLL Algorithm: Survey and Applications" edited by Phong
Q. Nguyen and Brigitte Vallée, Springer, Heidelberg, 2010.

See also Chapter 9 of M. R. Bremner, "Lattice Basis Reduction: An
Introduction to the LLL Algorithm and Its Applications" CRC Press, 2012.

# Examples
```jldoctest
julia> x = [0.3912641745333527; 0.5455179974014548; 0.1908698210882469];
julia> rationalapprox(x,1e4,Int64)
3-element Array{Rational{Int64},1}:
 43//110
  6//11
 21//110
```
"""
function rationalapprox(x::AbstractArray{<:Real,1},M,Ti=BigInt,verbose=false)
    # Follows the multivariate technique in Theorem 8 of the
    # Nguyen and Vallee paper above

    n = length(x)
    Q = (M/2^(1/4))^(1/n)
    C = Ti(ceil(Q^(n+1)))
    Cx = Ti.(round.(C*x))

    X = [one(Ti) zeros(Ti,1,n);
         [Cx     C*I] ]
    B,_ = lll(X)
    q = abs(B[1,1])
    maxErr = maximum(abs.(x-round.(x*q)/q))
    bound = sqrt(5)*2^(n/4 - 1)/Q /q
    if maxErr>bound
        warning("The max error in the approximation $(maxErr) is unexpectedly "*
                "greater than the bound $(bound).")
    end
    p = Ti.(round.(x*q));
    if verbose
        println("q        = $(q)\n"*
                "maxErr/q = $(maxErr)\n"*
                "bound/q  = $(bound)")
        display([x round.(x*q)/q abs.(x-round.(x*q)/q)])
    end
    return p//q
end
