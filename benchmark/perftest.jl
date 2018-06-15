if VERSION<=v"0.6.3"
    include("../src/LLLplus.jl")
else
    using Pkg
end

Pkg.add("PyPlot")
Pkg.add("BenchmarkTools")
Pkg.add("Docile")

using PyPlot
using BenchmarkTools
using LLLplus

include("lrtest.jl")

lrtest(5,2.^[2],[100],[Int32,Int64,Int128,Float64,BigInt,BigFloat],"rand")
savefig("benchmark/perfVsDataTypeN16.png") # run from root directory

lrtest(5,2.^[0:5;],[1],[Float64],"randn")
savefig("benchmark/perfVsNfloat32.png") # run from root directory
