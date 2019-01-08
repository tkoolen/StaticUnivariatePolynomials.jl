using StaticUnivariatePolynomials
using Test

import StaticUnivariatePolynomials: constant, derivative, integral

@testset "evaluation" begin
    @test Polynomial(1)(2) == 1
    @test Polynomial(2, 3, 4)(5) == 2 + 3 * 5 + 4 * 5^2
end

@testset "utility" begin
    p = Polynomial(2, 3, 4)
    @test constant(p) === 2
    @test zero(p) === Polynomial(0, 0, 0)
    @test transpose(p) === p
    @test conj(p) === p
    @test conj(Polynomial(1 + 2im, 3 + 4im, 5 - 6im)) === Polynomial(1 - 2im, 3 - 4im, 5 + 6im)
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

    p3 = Polynomial(ntuple(i -> 21 - i, Val(20)))
    derivative(p3)
    allocs = @allocated derivative(p3)
    @test allocs == 0
end

@testset "integral" begin
    p = Polynomial(ntuple(i -> 21 - i, Val(20)))
    c = 5
    P = integral(p, c)
    @test constant(P) == c
    P′ = derivative(P)
    for x in 0 : 0.1 : 1
        @test P′(x) ≈ p(x) atol=1e-14
    end
end
