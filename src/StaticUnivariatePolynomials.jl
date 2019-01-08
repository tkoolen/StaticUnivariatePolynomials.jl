module StaticUnivariatePolynomials

export
    Polynomial

struct Polynomial{N, T}
    coeffs::NTuple{N, T}
end

Polynomial(coeffs::Tuple) = Polynomial(promote(coeffs...))
Polynomial(coeffs...) = Polynomial(coeffs)

@inline function Polynomial{N}(p::Polynomial{M, T}) where {N, M, T}
    N < M && throw(InexactError(:Polynomial, Polynomial{N}, p))
    Polynomial((p.coeffs..., ntuple(_ -> zero(T), Val(N - M))...))
end

# Evaluation
(p::Polynomial{1})(x::Number) = p.coeffs[1] # evalpoly doesn't handle N = 1 case
@generated function (p::Polynomial{N})(x::Number) where N
    quote
        coeffs = p.coeffs
        @evalpoly(x, $((:(p.coeffs[$i]) for i = 1 : N)...))
    end
end

# Arithmetic
for op in [:+, :-]
    @eval begin
        function Base.$op(p1::Polynomial{N}, p2::Polynomial{N}) where N
            c1 = p1.coeffs
            c2 = p2.coeffs
            Polynomial(ntuple(i -> $op(c1[i], c2[i]), Val(N)))
        end
        function Base.$op(p1::Polynomial{N}, p2::Polynomial{M}) where {N, M}
            P = max(N, M)
            $op(Polynomial{P}(p1), Polynomial{P}(p2))
        end
    end
end

derivative(p::Polynomial{1}) = Polynomial(zero(p.coeffs[1]))
# @generated function derivative(p::Polynomial{N}) where N

# end

end # module
