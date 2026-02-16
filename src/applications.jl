"""
    integerfeasibility(A,d)

Given matrix `A` and vector `d`, return an integer vector `x` which satisfies the
system `Ax=d`.

    integerfeasibility(A,d,true)

If the third argument is present and is `true`, then as well as returning a
solution `x`, also return a matrix `xNull` of vectors in the null space of
`A` which could be added to the `x` vector to find a solution which
satisfies a constraint such as `0 .≤ x .≤ u`; see the paper below.

This is not a robust tool, just a demo.

"Solving A System Of Diophantine Equations With Bounds On The Variables" by
Karen Aardal, Cor Hurkens, and Arjen Lenstra in Integer Programming and
Combinatorial Optimization, 6th International IPCO Conference, vol 1412, pp
229-242, 1998. See http://softlib.rice.edu/pub/CRPC-TRs/reports/CRPC-TR98782.pdf

# Examples
```jldoctest
julia> A=[10 1 -9; 1 8 8]; xtrue=[0; 2; 9]; d=A*xtrue;

julia> integerfeasibility(A,d)
3-element Vector{Int64}:
 0
 2
 9

julia> A=[10 1.1 -9.1; 1 8 8]; d=A*xtrue;

julia> integerfeasibility(A,d)
3-element Vector{Float64}:
 0.0
 2.0
 9.0

julia> n=20;m=30; A = rand(-10:10,n,m); xtrue = rand(0:10,m); d=A*xtrue;

julia> sum(abs.(xtrue - integerfeasibility(A,d) ))
0

```
"""
function integerfeasibility(A::AbstractArray{Td,2},
                            d::AbstractArray{Td,1}, flag=false) where {Td<:Number}
    m,n=size(A)
    # N1=10^4
    # N2=10^5
    N1=10^4
    N2=10^5
    B = [[I           zeros(Td,n,1)];
         zeros(Td,1,n) N1;
         N2*A          -N2*d]

    Bhat,_ = lll(B)

    # the columns are not always sorted. A sort here helps sometimes.

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

This function tries first the technique in the `lagariasodlyzko` function,
and if it fails, a solution via `mdsubsetsum` is attempted.

It appears that this function can also solve some integer relations
problems. See the first example.

# Examples
```jldoctest
julia> a=[1.5;.5;0;.1;.2]; s=2.2; x,_=subsetsum(a,s,true); s-x'*a
A binary Lagarias-Odlyzko solution was found.
A solution was found via lagariasodlyzko
0.0

julia> a=[32771,65543,131101,262187,524387,1048759, # from Bremner p 117
          2097523,4195057,8390143,16780259,33560539,
          67121039,134242091,268484171,536968403];

julia> s=891221976; x,_=subsetsum(a,s,false); s-x'*a
0

julia> N=40;a=rand(1:2^BigInt(256),N);xtrue=rand(Bool,N); s=a'*xtrue;

julia> setprecision(BigFloat,300); x,_=subsetsum(a,s,false); s-x'*a
0

```
"""
function subsetsum(a::AbstractArray{Td,1},ss::Td,verbose=false) where {Td<:Number}

    xlo,binarylo = lagariasodlyzko(a,ss,0,true,verbose);
    if binarylo===true && ~any(ismissing.(xlo))
        verbose && print("A solution was found via lagariasodlyzko\n")
        return xlo,binarylo
    else
        xlo,binarylo = lagariasodlyzko(a,ss,0,false,verbose);
        if binarylo===true && ~any(ismissing.(xlo))
            verbose && print("A solution was found via lagariasodlyzko\n")
            return xlo,binarylo
        else
            x = mdsubsetsum(a,ss,.5,1)
            if ~ismissing(x)
                verbose && print("A solution was found via mdsubsetsum\n")
                return x,true
            elseif binarylo===false
                return xlo,binarylo
            end
        end
    end
    return missing
end


"""
    x = lagariasodlyzko(a,s)

For a vector of integers `a`, and an integer `s`, try to find a binary
vector `x` such that `x'*a=s`. We use the LLL algorithm to find the
solution. This is not a robust tool, just a demo.

This follows the technique described by Lagarias and Odlyzko  in 
"Solving Low-Density Subset Sum Problems"  in Journal of ACM, Jan 1985.
Code based on http://web.eecs.umich.edu/~cpeikert/lic15/lec05.pdf
We can likely get better results using techniques described and referenced in
https://www-almasty.lip6.fr/~joux/pages/papers/ToolBox.pdf

It's odd that permuting the `a` vector in the second example given below
causes the alg to often not find a binary solution. Apparently this is a
common oddity with lattice solvers.

# Examples
```jldoctest
julia> a=[1.5;.5;0;.1;.2]; s=2.2; x,_=lagariasodlyzko(a,s); s-x'*a
0.0

julia> a=[32771,65543,131101,262187,524387,1048759, # from Bremner p 117
          2097523,4195057,8390143,16780259,33560539,
          67121039,134242091,268484171,536968403];

julia> s=891221976; x,_=lagariasodlyzko(a,s); s-x'*a
0

julia> N=40;a=rand(1:2^BigInt(256),N);xtrue=rand(Bool,N); s=a'*xtrue; 

julia> setprecision(BigFloat,300); x,_=lagariasodlyzko(a,s); s-x'*a
0

```
"""
function lagariasodlyzko(a::AbstractArray{Td,1},ss::Td,b=0,enableFlip=true,
                         verbose=false) where {Td<:Number}
    # page numbers below refer to lecture note above

    n = length(a)

    Ti = getIntType(Td)
    if ss<sum(a)/2 && enableFlip
        flipflag=true
        s = sum(a)-ss
    else
        flipflag = false
        s = ss
    end
    if b≤0  # b is the "B" parameter defined at bottom of page 2
        nt = Ti(n)
        b = Ti(ceil(sqrt(nt*2^nt)))
    end
    # BB is the 'B' matrix defined also at bottom of page 2
    BB = [Matrix{Ti}(I,n,n+1); -b*a' b*s]
    B,_ = lll(BB)
    #B = l2(BB)

    # Pick off the interesting parts of the B matrix
    Bp = B[1:n,1:n+1]
    # Check if any of the columns satisfy the requirement
    sums = abs.(a'*Bp)
    ixMatch = findall(sums.≈s)
    if length(ixMatch)==0
        @printf("A match that is close according to `isapprox` was not found.\n")
        err,ixx = findmin(sums .-s)
        if abs(err)<1e-4
            @printf("A solution w err %e will be used. Check that it's correct \n",
                    Float64(err))
            ixMatch = findall(sums.-s .≈ err) # lazy hack to find ixMatch
        else
            xr = zeros(Ti,n+1,1).*missing
            return xr, false
        end
    end

    binarySolution = false
    xr = zeros(Ti,n+1,1)
    for ix=1:length(ixMatch)
        # Pick off a column of B
        x = Bp[:,ixMatch[ix][2]]
        if !(any(x.<0)|any(x.>1))
            
            binarySolution= true
            xr = flipflag ? mod.(x .+ 1,2) : x
        end
        x = -x
        if !(any(x.<0)|any(x.>1))
            binarySolution= true
            xr = flipflag ? mod.(x .+ 1,2) : x
        end
    end

    if binarySolution
        verbose && print("A binary Lagarias-Odlyzko solution was found.\n")
    else
        if maximum(a)<2^(n^2/2) && verbose
            density = n/maximum(log2.(a))
            @printf("The density (%4.2f) of 'a' is not as low as required ", density)
            @printf("(%4.2f) for Lagarias-Odlyzko to work. \n",2/n)
        end

        if length(ixMatch)>=1
            verbose && print("A non-binary Lagarias-Odlyzko solution was "*
                             "found; check that it's correct\n")
            xr=Bp[:,ixMatch[1][2]]
        else
            verbose && print("A Lagarias-Odlyzko solution was not found\n")
            xr= missing.*xr
        end
    end
    return xr, binarySolution
end


"""
    x = mdsubsetsum(a,sM,ratio=.5,Kpm=3)

For a vector of integers `a`, and an integer `sM`, try to find a binary
vector `x` such that `x'*a=sM` using the technique from "Multidimensional
subset sum problem" [1][2]. A major goal of the technique is to solve
problems in which there are about 50% ones in `x`; other ratios of ones to zeros
can be specified in `ratio`.  The thesis also suggests searching `Kpm=3`
values around the nominal k. This technique is related to that in
[`subsetsum`](@ref) in that both use the LLL algorithm.  This is not a
robust tool, just a demo.

[1] https://scholarworks.rit.edu/theses/64/

[2] https://pdfs.semanticscholar.org/21a7/c2f9ff29507f1153aefcca04d1cd308e45c0.pdf

# Examples
```jldoctest
julia> using Random, LLLplus

julia> a=[1.5;.5;0;.1;.2]; s=2.2; x=mdsubsetsum(a,s); s-x'*a
0.0

julia> a=[32771,65543,131101,262187,524387,1048759, # from Bremner p 117
          2097523,4195057,8390143,16780259,33560539,
          67121039,134242091,268484171,536968403];

julia> sM=891221976; x=mdsubsetsum(a,sM); sM-x'*a
0

julia> setprecision(BigFloat,300); Random.seed!(0);

julia> N=40;a=rand(1:2^BigInt(256),N);xtrue=rand(Bool,N); s=a'*xtrue;

julia> x=mdsubsetsum(a,s); s-x'*a
0
```
"""
function mdsubsetsum(a::AbstractArray{Td,1},sM::Td,ratio=.5,Kpm=3,
                     verbose=false) where {Td<:Number}

    Ti = getIntType(Td)
    n = Ti(length(a))

    # r is heuristic value for r from end of Section 11 of thesis above
    r = Ti(ceil(2^round(.65*log2(maximum(a)))))
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
        B,_ = lll(BB)
        #B = l2(BB)

        for nx = 1:n
            if abs(B[n+3,nx])==1 && B[n+2,nx]==0 && B[n+1,nx]==0 &&
                      all((B[1:n,nx] .==1) .| (B[1:n,nx] .==-1))
                xhat = abs.(B[1:n,nx] .- B[n+3,nx])/2
                return Ti.(xhat)
            end
        end
    end
    verbose && print("Solution not found")
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
example the Hanrot paper below.

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
3-element Vector{Rational{Int64}}:
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


"""
    spigotBBPvec(Td::Type{Tr},s,b,n,K) where {Tr<:Number}

This is an auxiliary function for `spigotBBP`. It calculates the numbers
that `spigotBBP` uses to check for a BBP spigot infinite series.

If you want to try to get BBP coefficients without the aid of the
`spigotBBP` function, below is code you can play with.

# Example
```jldoctest
julia> vv=LLLplus.spigotBBPvec(Float64,1,16,8,20);

julia> v=[vv;-pi];

julia> b=1_000_000;

julia> A=[I; b* v'];

julia> B,T=lll(A);

julia> B[1:end-2,1]'*vv -pi <1e-14 # The error is small
true
```
"""
function spigotBBPvec(Td::Type{Tr},s,b,n,K) where {Tr<:Number}
    v = Vector{Td}(undef,n)
    k = Td(0):K
    for j =1:n
        v[j]= sum(1 ./(b.^k .*(k.*n .+j).^s))
    end
    return v
end

"""
    spigotBBP(α::Td,s,b,n,K,verbose=false) where {Td}

Check for a BBP-style [1] infinite series for the constant `α`.  These are
"spigot" formulas that can be used to generate (for example) the millionth
digit of the constant `α` without learning the previous
digits. Specifically, given the constant `α`, and parameters `b`, `n`, and
`s`, look for a vector of numbers `a_1` through `a_n` that satisfies the
following equation:

``\\alpha= \\sum_{k=0}^\\infty \\frac{1}{b^k} \\left( \\frac{a_1}{(nk+1)^s} + \\ldots + \\frac{a_n}{(nk+n)^s} \\right)``

Because it's hard to sum to infinity, the sum is stopped at K.
If a formula is found, it is printed to the screen in LaTeX and the
coefficents `a` are returned as a vector.  An online LaTeX viewer
such as https://www.latex4technics.com/ may be helpful.

This is not a robust tool, just a demo. For example, there may be a 
problem with s≥2. See [2] for derivation of the technique used, and to 
check whether a formula you find is new.

[1] David Bailey, Peter Borwein, and Simon Plouffe. "On the rapid
computation of various polylogarithmic constants." Mathematics of
Computation 66.218 (1997): 903-913.
https://www.ams.org/journals/mcom/1997-66-218/S0025-5718-97-00856-9/

[2] David Bailey, "A Compendium of BBP-Type Formulas for Mathematical
Constants". https://www.davidhbailey.com//dhbpapers/bbp-formulas.pdf

# Example
```jldoctest
julia> spigotBBP(BigFloat(pi),1,16,8,45,true);
A solution was found w error -4.728672e-60. In LaTeX form it is
\\alpha= \\sum_{k=0}^\\infty \\frac{1}{16^k} \\left(\\frac{4}{8k+1}-\\frac{2}{8k+4}-\\frac{1}{8k+5}-\\frac{1}{8k+6}\\right)
```

Other examples without output:
```julia
spigotBBP(Float64(pi),1,-4,4,22,true);
spigotBBP(log(2),1,2,2,30,true);
spigotBBP(9*log(3),1,9,2,30,true);
spigotBBP(atan(2)*8,1,16,8,30,true);
spigotBBP(8*sqrt(2)*log(1+sqrt(2)),1,16,8,25,true);
```

There is a formula for pi^2 which the following command should find, but it
does not find it. In fact the technique doesn't seem to work at all for 
s>2; It's not obvious what the problem is
```julia
spigotBBP(BigFloat(pi)*pi,2,64,6,25,true);
```
"""
function spigotBBP(α::Td,s,b,n,K,verbose=false) where {Td}
    v = spigotBBPvec(Td,s,b,n,K)
    bN = 100_000_000
    av,_ = lagariasodlyzko(v,α,bN);
    if verbose && ~ismissing(av[1])
        @printf("A solution was found w error %e. In LaTeX form it is\n",av'*v-α)
        @printf("\\alpha= \\sum_{k=0}^\\infty \\frac{1}{%d^k} \\left(",b)
        for ix = 1:length(av)
            if abs(av[ix])>10*eps()
                if sign(av[ix])==1
                    ix>1 && @printf("+")
                else
                    @printf("-")
                end
                if s==1
                    @printf("\\frac{%d}{%dk+%d}",abs(Int64(av[ix])),n,ix)
                else
                    @printf("\\frac{%d}{(%dk+%d)^%d}",abs(Int64(av[ix])),n,ix,s)
                end
            end
        end
        @printf("\\right)\n")
    end
    if ismissing(av[1])
        verbose && @printf("A solution was not found.\n")
        return missing
    else
        return av
    end
end


"""
    minimalpolynomial(d,β::Td,p,verbose=false) where {Td<:Number}

Find the minimal polynomial for some element 'α' of a field given the degree
'd' of the polynomial, an approximation 'β' of 'α', and the number of
significant (decimal) digits of precision 'p' in 'β'. Returns two things, a
function, and a Bool indicating whether we think a true minimal polynomial
was found.  The ouptut function takes a single argument of the same type as
'β', and evaluates the estimated minimal polynomial over that argument.

Follows technique in the following paper: "A Gentle Tutorial for
Lattice-Based Cryptanalysis" by Joseph Surin and Shaanan Cohney,
https://eprint.iacr.org/2023/032

# Example
```jldoctest
julia> minimalpolynomial(3,8.70997594,8,true);
Inputs:
  d=3, β=8.709976, Td=Float64, p=8
Output:
  f(x) = x^3 - 21x^2 + 147x^1 + -348

julia> f,success = minimalpolynomial(3,8.70997594,9); success
false

```
"""
function minimalpolynomial(d,β::Td,p,verbose=false) where {Td<:Number}
    Ti = getIntType(Td)
    B = zeros(Ti,d+1,d+1)
    for ix = 1:d+1
        B[1,ix] = floor(10^p*β^(ix-1))
    end
    B[2:d+1,1:d] = Matrix{Ti}(I,d,d)
    Bp,_ = lll(B)
    a = Bp[2:d+1,1]

    # create f(x)
    function ff(x::Td)
        s = zero(Td)
        for ix = 1:d
            s += a[ix] * x^(ix-1)
        end
        s += x^d
        return s
    end

    # compare f(β) to β
    f = ff(β)
    success = true
    if abs(β/f) < 10^p
        success = false
    end

    if verbose
        @printf("Inputs:\n  d=%d, β=%f, Td=%s, p=%d\n",d,β,Td,p)
        @printf("Output:\n")
        @printf("  f(x) = x^%d",d)
        for ix = d:-1:2
            @printf(" %c %dx^%d",a[ix]<0 ? '-' : '+', abs(a[ix]),ix-1)
        end
        @printf(" + %d\n",a[1])
        if ~success
            @printf("The results look off\n")
            @printf("  For a high-precision β, we expect f(β)≈0, here f(β) = %f\n",f)
            @printf("  We expected to see abs(β/f(β)) > 10^%d, but instead saw ",p)
            @printf("abs(β/f(β)) = %f\n",abs(β/f))
            @printf("  We expected abs(β/f)/10^p>1, here abs(β/f)/10^%d = %f\n",
                    p,abs(β/f)/10^p)
            @printf("  Maybe β did not have %d digits of precision.\n",p)
            @printf("  or the minimal polynomial does not have degree %d.\n",d)
            x = abs.(Bp[:,2:end])[:]
            aveB = sum(x)/length(x)
            ratio = Bp[1,1]/aveB
            @printf("  When succesful, we expect the ratio between the first element of\n")
            @printf("  the shortest vector and the average of elements in other vectors\n")
            @printf("  to be small or less than 1. Here the ratio was %f \n",ratio)
        end
    end
    return ff,success
end


"""
    x = partitionintwo(a)

For a vector of integers `a`, check if the set `a` can be separated into two
sets that each have the same sum. I.e. try to find a binary vector `x` such
that for `y = .~ Bool.(x)` then `x'*a= y'*a. If no such partition is found,
`missing` is returned, otherwise `x` and `y` are returned.

# Examples
```jldoctest
julia> a = [3,1,1,2,2,1]; # from wikipedia page for Partition Problem

julia> x,y=partitionintwo(a); [x y]
6×2 BitMatrix:
 1  0
 0  1
 0  1
 0  1
 1  0
 0  1

```
"""
function partitionintwo(a::AbstractArray{Td,1}) where {Td<:Number}
    s = Td(sum(a)/2)

    x,_ = subsetsum(a,s)
    x = Bool.(x)
    if ismissing(x[1])
        return missing
    else
        y = .~ Bool.(x)
        sx = x'*a
        sy = y'*a

        if sx==sy
            return x,y
        else
            return missing
        end
    end
end
