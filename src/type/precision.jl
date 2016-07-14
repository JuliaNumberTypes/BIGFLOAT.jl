
function precision(x::BigFloat)  # precision of an object of type BigFloat
    return ccall((:mpfr_get_prec, :libmpfr), Clong, (Ptr{BigFloat},), &x)
end

precision(::Type{BigFloat}) = BIGFLOAT_PRECISION[end]  # precision of the type BigFloat itself

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
function setprecision{T}(f::Function, ::Type{T}, prec::Integer)
    old_prec = precision(T)
    setprecision(T, prec)
    try
        return f()
    finally
        setprecision(T, old_prec)
    end
end

setprecision(f::Function, precision::Integer) = setprecision(f, BigFloat, precision)
