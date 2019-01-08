using StaticUnivariatePolynomials
using Test

import StaticUnivariatePolynomials: derivative

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
    p3 = Polynomial(ntuple(i -> 21 - i, Val(20)))
    @test (p1 + p2)(7) == p1(7) + p2(7)
    allocs = @allocated p1 + p2
    @test allocs == 0
    @test (p1 + p3)(7) == p1(7) + p3(7)
    allocs = @allocated p1 + p3
    @test allocs == 0
end

@testset "derivative" begin
    p1 = Polynomial(3)
    @test derivative(p1) === Polynomial(0)

    p2 = Polynomial(1, 2, 3, 4)
    @test derivative(p2) === Polynomial(2, 2 * 3, 3 * 4)
end
