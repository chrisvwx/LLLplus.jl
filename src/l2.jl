lll(H::AbstractArray{Rational{Td},2},δ=.75,η=.51) where {Td<:Integer} =
    l2(H,Rational{Td},δ,η)

"""
    B = l2(H::AbstractArray{Td,2},TG::Type{Tg},δ=.75,η=.51) where {Td<:Number,Tg<:Number}

Do L2 lattice reduction of matrix `H` using optional parameters `δ` and `η`,
using a Gram matrix of type TG. The output is `B`, an LLL-reduced basis such
that `H*T = B`, where `T` is a unimodular matrix which is not returned.

Follows L2 from Damien Stehle, "Floating-Point LLL: Theoretical and
Practical Aspects", chapter 5 of "The LLL Algorithm, survey and
applications", Springer-Verlag, 2009.

The `l2` function is still in a raw (alpha) state and may change or be removed.
The `lll` function in this library also does LLL reduction, and functions on
a wider set of input types, such as complex data. When `l2` works, it is
generally faster than `lll` on small bases, say of dimensions less than 80. 

# Examples
```jldoctest
julia> H= [1 2; 3 4];B,_ = LLLplus.l2(H); B
┌ Warning: l2 is in a raw (alpha) state and may change. See the help text.
└ @ LLLplus src/l2.jl:45
2×2 Matrix{Int64}:
 1  -1
 1   1

julia> H= [.5 2; 3 4]; B,_= LLLplus.l2(H); B
┌ Warning: l2 is in a raw (alpha) state and may change. See the help text.
└ @ LLLplus src/l2.jl:45
2×2 Matrix{Float64}:
 1.5  -1.0
 1.0   2.0

julia> N=30;H = randn(N,N); B,T = LLLplus.l2(H);
┌ Warning: l2 is in a raw (alpha) state and may change. See the help text.
└ @ LLLplus src/l2.jl:45

```
"""
function l2(H::AbstractArray{Td,2},TG::Type{Tg}=Td,δ=.75,η=.51) where
    {Td<:Any,Tg<:Any}

    @warn "l2 is in a raw (alpha) state and may change. See the help text." maxlog=1 _file="src/l2.jl"
    
    if !(0.25 < δ < 1.0)
        error("δ must be between 1/4 and 1.");
    end
    if !(0.5 < η < sqrt(δ))
        error("η must be between 1/2 and sqrt(δ).");
    end
    n,d = size(H)
    Ti= LLLplus.getIntType(Td)
    Tf = float(Td)
    
    # Follows Stehle
    ϵ = .001  # ϵ = eps(Td) may be too small
    C = ϵ
    if Tf <: Complex
        @error "`l2` does not handle complex data; try `lll`."
        return
    end
    Tfe = real(Tf)
    l = precision(Tfe)
    left = d^2*( ((1+η)^2+ϵ)/(δ-ϵ^2) )^d * 2^(-l+10+C*d)
    right = minimum([ϵ η-.5 1-δ])
    # if left>right
    #     @warn "Precision of l=$l bits for $Tf is not sufficient "*
    #           "to satisfy Stehle's requirement [survey book, ch 5, thm 2]."
    # end
    if l<1.6*d
        @warn "Precision of l=$l bits does not even satisfy Stele's "*
        "asymptotic requirement of l≥1.6d." 
    end

    δ= Tfe(δ);
    η= Tfe(η)
    B = copy!(similar(H), H)
    G = Tg.(B)'*B
    ηb= (η+.5)/2
    δb= (δ+1)/2

    μ = Matrix{Tf}(I,n,d)
    r = zeros(Tf,n,d) 
    s = Vector{Tf}(undef,d)
    X = Vector{Ti}(undef,d)

    r[1,1] = G[1,1]

    κ=2;
    while κ<=d
        lazysizereduction!(ηb,κ,B,G,r,μ,s,X,n,d,Tg)

        κp=κ
        while κ>=2 && δb*real(r[κ-1,κ-1]) > real(s[κ-1])
            κ-=1
        end
        for i=1:κ-1
            μ[κ,i] = μ[κp,i]
            r[κ,i] = r[κp,i]
        end
        r[κ,κ] = s[κ]

        if κ!=κp
            bκp = B[:,κp]
            copyto!(B,n*κ+1,B,n*(κ-1)+1,n*length(κ:κp-1))
            B[:,κ] .= bκp
            
            if κ<d
                for i=κ:κp
                    @inbounds for d1=1:d
                        G[d1,i] = B[1,d1]*B[1,i]
                        for d2=2:d
                            G[d1,i] += B[d2,d1]*B[d2,i]
                        end
                    end
                end
                # maybe make G a `Hermitian` view of an underlying matrix,
                # set values using G.data[...] above, and skip the loop
                # below? see also discourse: https://bit.ly/36QguCG
                
                for i=κ:κp
                    for dx=1:d
                        G[i,dx] = G[dx,i]
                    end
                end
            end
        end
        κ+=1
    end
    T = Ti.(round.(inv(H)*B))
    return B,T
end

function lazysizereduction!(ηb,κ,B,G,r,μ,s,X,n,d,Tg)

    @label startCholesky
    # update the cholesky factorization
    i = κ
    for j = 1:i-1
        r[i,j] = G[i,j]
        for k = 1:j-1
            r[i,j] -= μ[j,k]*r[i,k]
        end
        μ[i,j] = r[i,j]/r[j,j]
    end
    s[1] = G[i,i]
    for j = 2:i
        s[j] = s[j-1]-μ[i,j-1]*r[i,j-1]
    end
    r[i,i] = s[i]

    μκ = view(μ,κ,1:κ-1)
    if maximum(abs.(μκ))≤ηb
        return
    else
        @inbounds for i = κ-1:-1:1
            X[i] = round(μ[κ,i])
            for j=1:i-1
                μ[κ,j]-=X[i]*μ[i,j]
            end
        end
        @inbounds for i = κ-1:-1:1
            for nx=1:n
                B[nx,κ] -=X[i]*B[nx,i]
            end
        end
        @inbounds for d1=1:d
            G[d1,κ] = B[1,d1]*B[1,κ]
            for d2=2:d
                G[d1,κ] += B[d2,d1]*B[d2,κ]
            end
        end
         for dx=1:d
            G[κ,dx] = G[dx,κ]
        end
    end
    @goto startCholesky
end
#=
"""
    B = l2turbo(H::AbstractArray{Td,2},TG::Type{Tg},δ=.75,η=.51) where
                         {Td<:Number,Tg<:Number}

A version of the `l2` function with a few calls to the `@turbo` macro
from the `LoopVectorization.jl` package.  See the `l2` help text

# Examples
```julia
julia> using LLLplus
julia> using LoopVectorization
julia> H= [1 2; 3 4];B = l2turbo(H)
┌ Warning: l2turbo is in a raw (alpha) state and may change. See the help text.
└ @ LLLplus ~/shared/LLLplus/src/l2.jl:42
2×2 Matrix{Int64}:
 1  -1
 1   1
```
"""
function l2turbo(H::AbstractArray{Td,2},TG::Type{Tg}=Td,δ=.75,η=.51) where
    {Td<:Number,Tg<:Number}

    @warn "l2turbo is in a raw (alpha) state and may change. See the help text." maxlog=1
    
    if !(0.25 < δ < 1.0)
        error("δ must be between 1/4 and 1.");
    end
    if !(0.5 < η < sqrt(δ))
        error("η must be between 1/2 and sqrt(δ).");
    end
    n,d = size(H)
    Ti= LLLplus.getIntType(Td)

    Tf = float(Td)
    
    # Follows Stehle 
    ϵ = .001  # ϵ = eps(Td) may be too small
    C = ϵ
    if Tf <: Complex
        @error "`l2turbo` does not handle complex data; try `lll`."
        return
    end
    Tfe = real(Tf)
    l = precision(Tfe)
    left = d^2*( ((1+η)^2+ϵ)/(δ-ϵ^2) )^d * 2^(-l+10+C*d)
    right = minimum([ϵ η-.5 1-δ])
    # if left>right
    #     @warn "Precision of l=$l bits for $Tf is not sufficient "*
    #           "to satisfy Stehle's requirement [survey book, ch 5, thm 2]."
    # end
    if l<1.6*d
        @warn "Precision of l=$l bits does not even satisfy Stele's "*
              "asymptotic requirement of l≥1.6d."
    end

    δ= Tfe(δ);
    η= Tfe(η)
    B = copy!(similar(H), H)
    G = Tg.(B)'*B
    ηb= (η+.5)/2
    δb= (δ+1)/2

    μ = Matrix{Tf}(I,n,d)
    r = zeros(Tf,n,d) 
    s = Vector{Tf}(undef,d)
    X = Vector{Ti}(undef,d)

    r[1,1] = G[1,1]

    κ=2;
    while κ<=d
        lazysizereductionAVX!(ηb,κ,B,G,r,μ,s,X,n,d,Tg)

        κp=κ
        while κ>=2 && δb*real(r[κ-1,κ-1]) > real(s[κ-1])
            κ-=1
        end
        for i=1:κ-1
            μ[κ,i] = μ[κp,i]
            r[κ,i] = r[κp,i]
        end
        r[κ,κ] = s[κ]

        if κ!=κp
            bκp = B[:,κp]
            copyto!(B,n*κ+1,B,n*(κ-1)+1,n*length(κ:κp-1))
            B[:,κ] .= bκp
            
            if κ<d
                @turbo for i=κ:κp
                    for d1=1:d
                        G[d1,i] = B[1,d1]*B[1,i]
                        for d2=2:d
                            G[d1,i] += B[d2,d1]*B[d2,i]
                        end
                    end
                end
                # maybe make G a `Hermitian` view of an underlying matrix,
                # set values using G.data[...] above, and skip the loop
                # below? see also discourse: https://bit.ly/36QguCG
                
                @inbounds for i=κ:κp
                    for dx=1:d
                        G[i,dx] = G[dx,i]
                    end
                end
            end
        end
        κ+=1
    end
    return B
end

function lazysizereductionAVX!(ηb,κ,B,G,r,μ,s,X,n,d,Tg)

    @label startCholesky
    # update the cholesky factorization
    i = κ
    @inbounds for j = 1:i-1
        r[i,j] = G[i,j]
        for k = 1:j-1
            r[i,j] -= μ[j,k]*r[i,k]
        end
        μ[i,j] = r[i,j]/r[j,j]
    end
    s[1] = G[i,i]
    @inbounds for j = 2:i
        s[j] = s[j-1]-μ[i,j-1]*r[i,j-1]
    end
    r[i,i] = s[i]

    μκ = view(μ,κ,1:κ-1)
    if maximum(abs.(μκ))≤ηb
        return
    else
        @inbounds for i = κ-1:-1:1
            X[i] = round(μ[κ,i])
            for j=1:i-1
                μ[κ,j]-=X[i]*μ[i,j]
            end
        end
        @turbo for nx=1:n
            for i = 1:κ-1
                B[nx,κ] -=X[i]*B[nx,i]
            end
        end
        @turbo for d1=1:d
            G[d1,κ] = B[1,d1]*B[1,κ]
            for d2=2:d
                G[d1,κ] += B[d2,d1]*B[d2,κ]
            end
        end
        for dx=1:d
            G[κ,dx] = G[dx,κ]
        end
    end
    @goto startCholesky
end
=#



"""
    dataTypeForGram(nbits,d)

Given a bit depth `nbits`, dimension `d` of a integer square matrix, return
the integer datatype needed to represent the Gram matrix in full
precision. This is conservative for most matrices, yet in some corner cases,
the deeper bit depth in the datatype returned is needed.

"""
function dataTypeForGram(nbits,d)
    # nbits*2-1 bits required to represent output of each multiply
    # additional log2(d) bits required for sums. So for each element of the
    # Gram matrix we need `bitsRequired` bits.
    bitsRequired = Int64(ceil(nbits*2-1+log2(d)))
    return intTypeGivenBitsRequired(bitsRequired)
end

"""
    dataTypeForGram(B::AbstractArray{Td,2}) where {Td<:Number}

Given a matrix B of eltype Td, return a datatype for representing the Gram
matrix of B. When Td<:Integer, we calculate the Gram matrix as BigInt, then
use its maximum value to get an `nbits` that to calculate a Gram data type
from `dataTypeForGram(nbits,d)`. Yes, it may be overkill to calculate the
Gram to find out the bit precision needed to represent it.

"""
function dataTypeForGram(B::AbstractArray{Td,2}) where {Td<:Number}
    if Td<:AbstractFloat
        return Float64  # maybe return Td?
    elseif Td<:Integer
        G = big.(B)'*B
        bitsRequired = Int64(ceil(log2(maximum(abs.(G)))))
        return intTypeGivenBitsRequired(bitsRequired)
    #elseif Complex
        #?
    else
        return Td
    end
end

"""
    intTypeGivenBitsRequired(bitsRequired)

Given the number of `bitsRequired`, return the smallest Integer type that
can represent that many bits, as well as a sign bit.

"""
function intTypeGivenBitsRequired(bitsRequired)
#    intTypes = [Int8, Int16, Int32, Int64,Int128,Int256,Int512,Int1024]
    intTypes = [Int8, Int16, Int32, Int64,Int128]
    intBits = sizeof.(intTypes)*8 .-1 # 1 bit for sign. could add 1 for spare
    ix = findfirst(bitsRequired .< intBits)
    return !isnothing(ix) ? intTypes[ix] : BigInt
end
