

# The default precision in significand bits for the type family `BigFloat`,
# which is alterable; this known with certainty only when called.
precision(::Type{BigFloat}) = BIGFLOAT_PRECISION[end]

# precision in signficand bits of an object of type `BigFloat`
precision{P}(x::BigFloat{P}) = P

"""
    setprecision([T=BigFloat,] precision::Int)

Set the precision (in bits) to be used for `T` arithmetic.
"""
function setprecision(::Type{BigFloat}, precision::Int)
    if precision < 2
        throw(DomainError())
    end
    BIGFLOAT_PRECISION[end] = precision
end

setprecision(precision::Int) = setprecision(BigFloat, precision)


"""
    setprecision(f::Function, [T=BigFloat,] precision::Integer)

Change the `T` arithmetic precision (in bits) for the duration of `f`.
It is logically equivalent to:

    old = precision(BigFloat)
    setprecision(BigFloat, precision)
    f()
    setprecision(BigFloat, old)

Often used as `setprecision(T, precision) do ... end`
"""
function setprecision(f::Function, ::Type{BigFloat}, prec::Integer)
    old_prec = precision(BigFloat)
    setprecision(BigFloat, prec)
    try
        return f()
    finally
        setprecision(BigFloat, old_prec)
    end
end

setprecision(f::Function, precision::Integer) = setprecision(f, BigFloat, precision)
