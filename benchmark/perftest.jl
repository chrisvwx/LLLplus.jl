include("../src/LLLplus.jl")
Pkg.add("PyPlot")
Pkg.add("CPUTime")

using PyPlot
using LLLplus
using CPUTime

include("lrtest.jl")


lrtest(50,2.^[4;],[100;],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
lrtest(50,2.^[4;],[100;],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
savefig("benchmark/perfVsDataTypeN16.svg") # run from root directory

lrtest(50,2.^[0:6;],[1;],[Float64],"randn")
lrtest(50,2.^[0:6;],[1;],[Float64],"randn")
savefig("benchmark/perfVsNfloat32.svg") # run from root directory
