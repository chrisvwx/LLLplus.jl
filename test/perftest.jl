using Plots
using LLLplus
using Printf
using Random
using LinearAlgebra
using DoubleFloats
using Quadmath
gr();

include("lrtest.jl")

# getIntType is used to indicate what type we want the unimodular integer
# matrix to have. We need to add methods of getIntType to handle the
# DoubleFloats and Quadmath, since they are external packages that LLLplus
# doesn't know about. If this is type piracy, I think it's acceptable.
import LLLplus.getIntType
getIntType(Td::Type{Tr}) where {Tr<:Float128} = Int128
getIntType(Td::Type{Tr}) where {Tr<:Double64} = Int128

lrtest(100,2 .^[4],[100],[Int32,Int64,Int128,Float32,Float64,Double64,Float128,BigInt,BigFloat],"rand")
savefig("perfVsDataType.png")
#display(current())

# lrtest(40,2 .^[1:8;],[1],[Float64],"randn")
# savefig("perfVsNfloat64.png")

