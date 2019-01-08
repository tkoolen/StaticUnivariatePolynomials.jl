using StaticUnivariatePolynomials
using Test

@testset "evaluation" begin
    @test Polynomial(1)(2) == 1
    @test Polynomial(2, 3, 4)(5) == 2 + 3 * 5 + 4 * 5^2
end

@testset "padding" begin
    p = Polynomial(2, 3, 4)
    q = Polynomial{20}(p)
    @test p(5) == q(5)
    allocs = @allocated Polynomial{20}(p)
    @test allocs == 0
end

@testset "arithmetic" begin
    p1 = Polynomial(2, 3, 4)
    p2 = Polynomial(6, 7, 8)
    @test (p1 + p2)(7) == p1(7) + p2(7)
    allocs = @allocated p1 + p2
    @test allocs == 0

    p3 = Polynomial(ntuple(_ -> rand(), Val(20)))
    @test (p1 + p3)(7) == p1(7) + p3(7)
end
