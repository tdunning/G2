using G2
using DataStructures

using Test

@testset "entropy" begin
    @test G2.entropy([1,1]) ≈ log(2)
    @test G2.entropy([1,1,1,1]) ≈ log(4)
end

@testset "g2" begin
    @test G2.g2([1 1;1 1]) == 0
    @test G2.g2([1 0;0 1]) ≈ 4 * log(2)
    @test G2.g2([10 0;0 10]) ≈ 40 * log(2)
    @test G2.g2([10 1; 2 3]) ≈ 4.562613816026058
end

@testset "g2root" begin
    @test G2.g2root([1 0; 0 1]) ≈ 2 * sqrt(log(2))
    @test G2.g2root([0 1; 1 0]) ≈ -2 * sqrt(log(2))
end

@testset "compare" begin
    r = G2.compare(counter("abcaba"), counter("abcbccdefdefdefdef"))
    @test r['f'] ≈ -1.59922 atol=1e-5
    @test r['a'] ≈ 2.36327 atol=1e-5
    @test r['c'] ≈ 0.0
end
