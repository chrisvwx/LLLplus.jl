# --------------
# Initialization
# --------------
include("../src/LLLplus.jl")

using LLLplus
using Base.Test

# --------------
# tests with small matrices
# --------------
println("tests with small matrices\n")

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
@test sum(abs(B-B1))==0
@test sum(abs(T-T1))==0

# --------------
# test norm for matrix from http://www.latticechallenge.org/
# --------------
m=200
println("Testing now with $(m)x$(m) matrix from latticechallenge.org.")
println("All the column norms should be 30.")
mat = readdlm("challenge-$(m).mod",Int64)

nrms = zeros(m,1)
for ix=1:m
    nrms[ix] = norm(mat[:,ix])
end
println("max col-norm of input is $(maximum(nrms))")

@time (B,T) = lll(mat);
nrms = zeros(m,1)
for ix=1:m
    nrms[ix] = norm(B[:,ix])
end
maxlll = maximum(nrms)
println("max col-norm of lll-reduced basis is $(maxlll)")
@test maxlll<=30+1e-6

@time (B,T) = seysen(mat);
nrms = zeros(m,1)
for ix=1:m
    nrms[ix] = norm(B[:,ix])
end
maxSeysen = maximum(nrms)
println("max column norm of seysen-reduced basis is $(maxSeysen)")
@test maxSeysen<=30+1e-6

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

# Test sphere decoder
Ns = 100000;
N = 2;
println("\nTesting sphere decoder on $(Ns) samples of $(N)x$(N) BPSK system...")
NN = randn(N,Ns)/sqrt(100);
H = [.6 -2.6; -.5 .95];
C = [-1,1];
Z = rand(1:2,N,Ns);
X = C[Z];
Y = H*X+NN;
@time Zt = hard_sphere(Y,H,2);
@time Zt = hard_sphere(Y,H,2);
Zh = Zt+1;
errRate = sum(abs(Z-Zh))/Ns;
println("Error Rate is $(errRate). It should be zero or very small.\n")
