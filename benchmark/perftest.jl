using Plots
using LLLplus
using Printf
using Random
using LinearAlgebra
using DoubleFloats
using Quadmath

include("lrtest.jl")

lrtest(40,2 .^[7],[100],[Int32,Int64,Int128,Float32,Float64,Float128,Double64,BigInt,BigFloat],"rand")
savefig("perfVsDataType.png")

lrtest(40,2 .^[0:8;],[1],[Float64],"randn")
savefig("perfVsNfloat64.png")
