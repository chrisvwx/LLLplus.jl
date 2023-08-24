## LLLplus README

```@meta
CurrentModule = LLLplus
```

LLLplus provides lattice tools such as
Lenstra-Lenstra-Lovász (LLL) lattice reduction which are of practical and
theoretical use in cryptography, digital communication, integer
programming, and more.
This package is experimental and not a robust tool; use at your own
risk :-)

LLLplus has functions for [LLL](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm),
[Seysen](http://link.springer.com/article/10.1007%2FBF01202355), and
[Hermite-Korkine-Zolotarev](http://www.cas.mcmaster.ca/~qiao/publications/ZQW11.pdf)
lattice reduction
techniques. [Brun](https://archive.org/stream/skrifterutgitavv201chri#page/300/mode/2up)
integer relations is included in the form of lattice
reduction. Solvers for the [shortest
vector](https://en.wikipedia.org/wiki/Lattice_problem#Shortest_vector_problem_(SVP))
and the [closest
vector](https://en.wikipedia.org/wiki/Lattice_problem#Closest_vector_problem_.28CVP.29)
problems are also included; for more see the help text for the `lll`,
`seysen`, `hkz`, `brun`, `svp`, and `cvp` functions. Several toy (demo)
functions are also included; see the  `subsetsum`, `minimalpolynomial`,
`integerfeasibility`, `rationalapprox`, and  `spigotBBP` functions.


### Examples

Each function contains documentation and examples available via Julia's
built-in documentation system (try `?lll` or `@doc(lll)`). Documentation
for all functions is [available](https://christianpeel.github.io/LLLplus.jl/dev). A tutorial notebook is
found in the `docs` directory or on
[nbviewer](https://nbviewer.jupyter.org/github/christianpeel/LLLplus.jl/blob/master/docs/LLLplusTutorial.ipynb).

Here are a few examples of using the functions in the
package on random lattices.

```julia
Pkg.add("LLLplus")
using LLLplus

# do lattice reduction on a matrix with randn entries
N = 40;
H = randn(N,N);
Bbrun,_ = brun(H);
Blll,_ = lll(H);
Bseysen,_ = seysen(H);
Bhkz,_ = hkz(H);

# check out the CVP solver
Q,Rtmp=qr(H); R = UpperTriangular(Rtmp);
u=Int.(rand(0:1e10,N));
y=H*u+rand(N)/100;
uhat=cvp(Q'*y,R);
sum(abs.(u-uhat))
```

### Execution Time results

To give a flavor of the behavior of the functions in LLLplus,
we show execution time for several built-in datatypes (Int32,
Int64, Int128, Float32, Float64, BitInt, and BigFloat) as well as type
from external packages (Float128 from
[Quadmath.jl](https://github.com/JuliaMath/Quadmath.jl) and Double64
from [DoubleFloat.jl](https://github.com/JuliaMath/DoubleFloats.jl))
which are used to generate 100 16x16 matrices with elements uniformly
distributed over `-100` to `100`. The figure shows average execution
time when using these matrices as input lattice bases for several
functions from LLLplus. See `test/perftest.jl` for the code to
regenerate the figure and for another line of code that generates a
figure of execution time versus basis dimension.

![Time vs data type](assets/perfVsDataType.png)

### Notes

The 2020 [Simons Institute lattice](https://simons.berkeley.edu/programs/lattices2020)
workshop, a
[survey paper by Wuebben](http://www.ant.uni-bremen.de/sixcms/media.php/102/10740/SPM_2011_Wuebben.pdf), and the
[monograph by Bremner](https://www.amazon.com/Lattice-Basis-Reduction-Introduction-Applications/dp/1439807027) 
were helpful in writing the tools in LLLplus
and are good resources for further study. If you are trying to break
one of the [Lattice Challenge](http://www.latticechallenge.org)
records or are looking for robust, well-proven lattice tools, look at
[fplll](https://github.com/fplll/fplll). Also, for many
number-theoretic problems the
[Nemo.jl](https://github.com/Nemocas/Nemo.jl) package is appropriate;
it uses the [FLINT](http://flintlib.org/) C library to do LLL
reduction on Nemo-specific data types.  Finally, no number theorists
have worked on LLLplus; please treat the package as experimental.
