

# The default precision in significand bits for the type family `BigFP`,
# which is alterable; this known with certainty only when called.
precision(::Type{BigFP}) = BIGFLOAT_PRECISION[end]

# precision in signficand bits of an object of type `BigFP`
precision{P}(x::BigFP{P}) = P

"""
    setprecision([T=BigFP,] precision::Int)

Set the precision (in bits) to be used for `T` arithmetic.
"""
function setprecision(::Type{BigFP}, precision::Int)
    if precision < 2
        throw(DomainError())
    end
    BIGFLOAT_PRECISION[end] = precision
end

setprecision(precision::Int) = setprecision(BigFP, precision)


"""
    setprecision(f::Function, [T=BigFP,] precision::Integer)

Change the `T` arithmetic precision (in bits) for the duration of `f`.
It is logically equivalent to:

    old = precision(BigFP)
    setprecision(BigFP, precision)
    f()
    setprecision(BigFP, old)

Often used as `setprecision(T, precision) do ... end`
"""
function setprecision(f::Function, ::Type{BigFP}, prec::Integer)
    old_prec = precision(BigFP)
    setprecision(BigFP, prec)
    try
        return f()
    finally
        setprecision(BigFP, old_prec)
    end
end

setprecision(f::Function, precision::Integer) = setprecision(f, BigFP, precision)
