"""
Main module for `LLLplus.jl` -- lattice reduction and related tools for Julia.

As an example of the functions in the package, see [`lll`](@ref), which does
Lenstra–Lenstra–Lovász lattice reduction of a matrix.

"""
module LLLplus

using LinearAlgebra
using Printf

export
    lll,
    cvp,
    svp,
    brun,
    gauss,
    sizereduction,
    seysen,
    vblast,
    subsetsum,
    integerfeasibility,
    rationalapprox,
    hardsphere, hard_sphere,
    issizereduced,islllreduced,orthogonalitydefect,hermitefactor,seysencond

include("lll.jl")          # lll, gauss, sizereduction
include("cvp.jl")          # cvp, svp
include("brun.jl")
include("seysen.jl")
include("vblast.jl")
include("applications.jl") # subsetsum, integerfeasibility, rationalapprox
include("utilities.jl")

include("hard_sphere.jl")  # may be deprecated in future

roundf(r::Td) where {Td<:Complex} = round(real(r)) + im*round(imag(r));
roundf(r) = round(r);

end # LLLplus module
