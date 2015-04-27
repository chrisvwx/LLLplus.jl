# LLLplus Performance results

On this page we give a few performance results from tests run on
Travis-CI during normal CI tests. In the tests we time execution of the
lattice-reduction functions, average the results over multiple random
matrices, and show results as a function of the size of the matrix and
of the data type. 

We first show how the time varies with matrix size (1,2,4,...64); the
vertical axis shows execution time on a logarithmic scale; the x-axis
is also logarithmic. The generally linear nature of the LLL curve supports
the polynomial-time nature of the algorithm. Each data point
is the average of execution time of 50 runs of a lattice-reduction
technique, where the matrices used were generated using *randn* to
emulate unit-variance Gaussian-distributed values.
![Time vs matrix size](perfVsNfloat32.png)

The lattice-reduction techniques work for Intgeger (Int64),
FloatingPoint (Float64), BigInt, and BigFloat. The vertical axis in
the next figure is a logarithmic representation of execution time as
in the previous figure. In the horizontal axis, the values 1..6
represent Int32, Int64, Int128, Float64, BitInt, and BigFloat
datatypes which are used to generate 50 16x16 matrices, over which
execution time for the lattic reduction techniques is averaged.
![Time vs data type](perfVsDataTypeN16.png)

#### Possible updates to the tests

* Basic implementations of the
  [Orthogonality defect](https://en.wikipedia.org/wiki/Lattice_reduction)
  and Hermite factor are in *lrtest.jl* now; these metrics could be
  explicit outputs for certain scenarios.
* The [SVP](http://www.latticechallenge.org/svp-challenge/) Challenge
  and the
  [Ideal](http://www.latticechallenge.org/ideallattice-challenge/)
  Lattice challenge have code to generate lattices for the respective
  contests which could be used or duplicated to make challenging
  tests. The main [Lattice](http://www.latticechallenge.org/)
  Challenge also lists references which could be used to replicate
  tests.
* I'd like to make better-looking plots; Gadfly is one option.
  
