module LLLplus

if VERSION<=v"0.6.3"
    using Base.BLAS
#    using ArrayViews
else
    using LinearAlgebra
end

# if VERSION < v"0.4.0-dev"
#     using Docile
# end
# @docstrings

export
    lll,
    seysen,
    vblast,
    hard_sphere

include("lll.jl")
include("seysen.jl")
include("vblast.jl")
include("hard_sphere.jl")

end # LLLplus module
