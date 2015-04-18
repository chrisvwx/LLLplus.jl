#
# Include LLLplus code  --------------
#
include("../src/LLLplus.jl")

using LLLplus

#
# Run tests  --------------
#

println("In all the following tests, the first time includes the "*
        "JIT compilation,\n while for the second execution the "*
        "compilation is done and the time should be faster.")
println(" ")


# Time LLL decomposition of a 1000x1000 real matrix with randn entries 
N = 1000;
H = randn(N,N);
println("Testing LLL on $(N)x$(N) real matrix...")
@time (B,T,Q,R) = lll(H);
@time (B,T,Q,R) = lll(H);

# Time LLL, Seysen, VBLAST decompositions of a 10x10 complex matrix with
# randn entries
N = 10;
H = randn(N,N) + im*randn(N,N);
println("Testing LLL on $(N)x$(N) complex matrix...")
@time (B,T,Q,R) = lll(H);
@time (B,T,Q,R) = lll(H);
println("Testing Seysen on $(N)x$(N) complex matrix...")
@time (T,H_red,H_red_dual) = seysen(H);
@time (T,H_red,H_red_dual) = seysen(H);
println("Testing VBLAST on $(N)x$(N) complex matrix...")
@time (W,P,B) = vblast(H);
@time (W,P,B) = vblast(H);

# Test sphere decoder
Ns = 100000;
N = 2;
println("Testing sphere decoder on $(Ns) samples of $(N)x$(N) BPSK system...")
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
println("Error Rate is $(errRate). It should be zero or very small.")
