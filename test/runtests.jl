# --------------
# Initialization
# --------------
using LLLplus
using Test
using LinearAlgebra
using Documenter

# --------------
# tests with small matrices
# --------------
@testset "tests with small matrices..." begin

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
end

# --------------
# Run tests with random matrices
# --------------

@testset "random matrices" begin
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
end

# Test CVP
@testset "CVP" begin
    N=500;
    println("\nTesting cvp on $(N)x$(N) BPSK system...")
    H=randn(N,N);
    Q,Rtmp=qr(H); R=UpperTriangular(Rtmp);
    u=Int.(rand([-1,1],N));
    y=H*u+rand(N)/100;
    uhat=cvp(Q'*y,R,Val(false),-1,1);
    errRate = sum(abs.(u-uhat))
    println("Error Rate is $(errRate). It should be zero or very small.\n")
end

# --------------
# doctests
# --------------
if VERSION >= v"1.1.0"
    println("\nRunning doctests...")
    DocMeta.setdocmeta!(LLLplus, :DocTestSetup,
                        :(using LLLplus, LinearAlgebra,Random); recursive=true)
    doctest(LLLplus)
end
