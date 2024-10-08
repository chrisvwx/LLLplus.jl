using StaticArrays, Test

@testset "Compatibility with StaticArrays/immutable arrays" begin
    H = [0.9628813734720784 0.23481870903463398 0.09437229281315995; 
         0.6398958697159318 0.19228011953230373 0.7960518809701681; 
         0.9766849678971312 0.6302875815478448 0.4788314270985644]
    sH = SMatrix{3, 3, Float64, 9}(H)
    @test lll(H)[1] ≈ lll(sH)[1]
    @test l2(H)[1] ≈ l2(sH)[1]
    @test brun(H)[1] ≈ brun(sH)[1]
    @test seysen(H)[1] ≈ seysen(sH)[1]
end