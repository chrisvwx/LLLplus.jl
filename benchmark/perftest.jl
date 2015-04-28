include("../src/LLLplus.jl")
Pkg.add("PyPlot")
Pkg.add("CPUTime")

using PyPlot
using LLLplus
using CPUTime

include("lrtest.jl")


lrtest(100,2.^[4;],[100;],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
lrtest(100,2.^[4;],[100;],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
savefig("benchmark/perfVsDataTypeN16.png") # run from root directory

lrtest(100,2.^[0:6;],[1;],[Float64],"randn")
lrtest(100,2.^[0:6;],[1;],[Float64],"randn")
savefig("benchmark/perfVsNfloat32.png") # run from root directory
