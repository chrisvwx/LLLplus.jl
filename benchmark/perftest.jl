include("../src/LLLplus.jl")

using PyPlot
using LLLplus
using USERTime

include("lrTest.jl")


lrTest(50,2.^[4;],[100;],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
lrTest(50,2.^[4;],[100;],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
savefig("benchmark/perfVsDataTypeN16.png") # run from root directory

lrTest(50,2.^[0:6;],[100;],[Float64],"rand")
lrTest(50,2.^[0:6;],[100;],[Float64],"rand")
savefig("benchmark/perfVsNfloat32.png") # run from root directory
