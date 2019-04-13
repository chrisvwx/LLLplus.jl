# to test code, run following in "test" directory
# julia> include("runtests.jl")

# --------------
# Initialization
# --------------
using LLLplus
using Test
using LinearAlgebra
using DelimitedFiles

# --------------
# tests with small matrices
# --------------
println("tests with small matrices...")

# Matrix from http://home.ie.cuhk.edu.hk/~wkshum/wordpress/?p=442
H =[1   9   1   2;
    1   8   8   3;
    7   4   5   1;
    2   6   7   1];
B1=[2   3  -2  -4
    3  -1   2   1
    1   1   6  -4
    1   3  -1   3];
T1=[0   0   1  -1
    0   1  -1   0
    0   0   0   1
    1  -3   3  -2];
(B,T) = lll(H);
@test sum(abs.(B-B1))==0
@test sum(abs.(T-T1))==0
println("...done")


# --------------
# Run tests with random matrices
# --------------

println("\nIn all the following tests, the first time includes the "*
        "JIT compilation; \nfor the second execution the "*
        "compilation is already done and the time\n"*
        "should be faster.\n")

# Time LLL decomposition of a 1000x1000 real matrix with randn entries
N = 1000;
H = randn(N,N);
println("Testing LLL on $(N)x$(N) real matrix...")
@time (B,T,Q,R) = lll(H);
@time (B,T,Q,R) = lll(H);

# Time LLL, Seysen, VBLAST decompositions of a 10x10 complex matrix with
# randn entries
N = 10;
println("\nTesting LLL on $(N)x$(N) complex matrix...")
H = randn(N,N) + im*randn(N,N);
@time (B,T,Q,R) = lll(H);
@time (B,T,Q,R) = lll(H);
println("Testing Seysen on same $(N)x$(N) complex matrix...")
@time (B,T) = seysen(H);
@time (B,T) = seysen(H);
println("Testing VBLAST on same $(N)x$(N) complex matrix...")
@time (W,P,B) = vblast(H);
@time (W,P,B) = vblast(H);
println("Testing Brun on real part of same $(N)x$(N) matrix...")
@time (B,T) = brun(real.(H));
@time (B,T) = brun(real.(H));

# Test sphere decoder
Ns = 100000;
N = 3;
println("\nTesting sphere decoder on $(Ns) samples of $(N)x$(N) BPSK system...")
NN = randn(N,Ns)/sqrt(100);
H = randn(N,N);
C = [-1,1];
Z = rand(1:2,N,Ns);
X = C[Z];
Y = H*X+NN;
@time Xt = hardsphere(Y,H,2);
@time Xt = hardsphere(Y,H,2);
errRate = sum(abs.(X-Xt))/Ns;
println("Error Rate is $(errRate). It should be zero or very small.\n")

# --------------
# test norm for matrix from http://www.latticechallenge.org/
# --------------

println("Testing now with 200x200 matrix from latticechallenge.org.")
println("The min norm of the input should be 30, min norm of the "*
        "reduced bases should be smaller.")
mat = readdlm("challenge-200.mod",Int64) #run from parent directory
mat = Matrix(mat');
N = size(mat,2)

nrms = zeros(N,1);
for ix=1:N
    nrms[ix] = norm(mat[:,ix])
end
println("min norm of input is $(minimum(nrms))")

@time (B,T) = lll(mat);
nrms = zeros(N,1);
for ix=1:N
    nrms[ix] = norm(B[:,ix]);
end
mlll=minimum(nrms)
println("min norm of lll-reduced basis is $(mlll)")
@test mlll<=30+1e-6

@time (B,T) = seysen(mat);
nrms = zeros(N,1)
for ix=1:N
    nrms[ix] = norm(B[:,ix])
end
mSeysen = minimum(nrms)
println("min norm of seysen-reduced basis is $(mSeysen)")
@test mSeysen<=30+1e-6

