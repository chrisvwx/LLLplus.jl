using PyPlot
using BenchmarkTools
using LLLplus
using Printf
using Random
using LinearAlgebra

include("lrtest.jl")

lrtest(5,2 .^[2],[100],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
savefig("benchmark/perfVsDataTypeN16.png") # run from root directory

lrtest(5,2 .^[0:5;],[1],[Float64],"randn")
savefig("benchmark/perfVsNfloat32.png") # run from root directory
