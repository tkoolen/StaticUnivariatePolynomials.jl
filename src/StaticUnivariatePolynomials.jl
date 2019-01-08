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
    coeffs = (p.coeffs..., ntuple(_ -> zero(T), Val(N - M))...)
    Polynomial(coeffs)
end

# Utility
constant(p::Polynomial) = p.coeffs[1]
Base.zero(::Type{Polynomial{N, T}}) where {N, T} = Polynomial(ntuple(_ -> zero(T), Val(N)))
Base.zero(p::Polynomial) = zero(typeof(p))
Base.conj(p::Polynomial) = Polynomial(map(conj, p.coeffs))
Base.transpose(p::Polynomial) = p

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
        # Two Polynomials
        function Base.$op(p1::Polynomial{N}, p2::Polynomial{N}) where N
            c1 = p1.coeffs
            c2 = p2.coeffs
            ntuple(Val(N)) do i
                $op(c1[i], c2[i])
            end |> Polynomial
        end
        function Base.$op(p1::Polynomial{N}, p2::Polynomial{M}) where {N, M}
            P = max(N, M)
            $op(Polynomial{P}(p1), Polynomial{P}(p2))
        end

        # Polynomial and Number
        Base.$op(p::Polynomial, c::Number) = Polynomial($op(constant(p), c), Base.tail(p.coeffs)...)
        Base.$op(c::Number, p::Polynomial) = Polynomial($op(c, constant(p)), map($op, Base.tail(p.coeffs))...)

        # Unary ops
        Base.$op(p::Polynomial) = Polynomial(map($op, p.coeffs))
    end
end

# Calculus
derivative(p::Polynomial{1}) = Polynomial(zero(p.coeffs[1]))
function derivative(p::Polynomial{N}) where N
    ntuple(Val(N - 1)) do i
        i * p.coeffs[i + 1]
    end |> Polynomial
end

function integral(p::Polynomial{N}, c) where N
    tail = ntuple(Val(N)) do i
        p.coeffs[i] / i
    end
    T = eltype(tail)
    Polynomial((T(c), tail...))
end

end # module
