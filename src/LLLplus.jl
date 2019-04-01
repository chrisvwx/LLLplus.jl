module LLLplus

using LinearAlgebra
using Printf

export
    lll,
    cvp,
    svp,
    gauss,
    seysen,
    vblast,
    subsetsum,
    integerfeasibility,
    hard_sphere

include("lll.jl")          # lll, gauss
include("cvp.jl")          # cvp, svp
include("seysen.jl")
include("vblast.jl")
include("applications.jl") # subsetsum, integerFeasibility

include("hard_sphere.jl")  # may be deprecated in future

roundf(r::Td) where {Td<:Complex} = round(real(r)) + im*round(imag(r));
roundf(r) = round(r);

end # LLLplus module
