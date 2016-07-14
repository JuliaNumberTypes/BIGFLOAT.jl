
widen(::Type{Float64}) = BigFloat
widen(::Type{BigFloat}) = BigFloat
widen{P}(::Type{BigFloat{P}}) = BigFloat{2*P+24}


maxintfloat(::Type{BigFloat}) = BigFloat(2)^precision(BigFloat)
maxintfloat{P}(x::BigFloat{P}) = BigFloat{P}(2)^precision(x)

signbit{P}(x::BigFloat) = ccall((:mpfr_signbit, :libmpfr), Int32, (Ptr{BigFloat{P}},), &x) != 0

function copysign{P}(x::BigFloat{P}, y::BigFloat{P})
    z = BigFloat{P}()
    ccall((:mpfr_copysign, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end

function exponent{P}(x::BigFloat{P})
    if x == 0 || !isfinite(x)
        throw(DomainError())
    end
    # The '- 1' is to make it work as Base.exponent
    return ccall((:mpfr_get_exp, :libmpfr), Clong, (Ptr{BigFloat{P}},), &x) - 1
end

function frexp{P}(x::BigFloat{P})
    z = BigFloat{P}()
    c = Clong[0]
    ccall((:mpfr_frexp, :libmpfr), Int32, (Ptr{Clong}, Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Cint), c, &z, &x, ROUNDING_MODE[end])
    return (z, c[1])
end

function significand{P}(x::BigFloat{P})
    z = BigFloat{P}()
    c = Clong[0]
    ccall((:mpfr_frexp, :libmpfr), Int32, (Ptr{Clong}, Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Cint), c, &z, &x, ROUNDING_MODE[end])
    # Double the significand to make it work as Base.significand
    ccall((:mpfr_mul_si, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Clong, Int32), &z, &z, 2, ROUNDING_MODE[end])
    return z
end


function isinteger{P}(x::BigFloat{P})
    return ccall((:mpfr_integer_p, :libmpfr), Int32, (Ptr{BigFloat{P}},), &x) != 0
end


function isinf{P}(x::BigFloat{P})
    return ccall((:mpfr_inf_p, :libmpfr), Int32, (Ptr{BigFloat{P}},), &x) != 0
end

function isnan{P}(x::BigFloat{P})
    return ccall((:mpfr_nan_p, :libmpfr), Int32, (Ptr{BigFloat{P}},), &x) != 0
end

isfinite{P}(x::BigFloat{P}) = !isinf(x) && !isnan(x)

@eval typemax(::Type{BigFloat}) = $(BigFloat( Inf))
@eval typemin(::Type{BigFloat}) = $(BigFloat(-Inf))
@eval typemax{P}(::Type{BigFloat{P}}) = $(BigFloat{P}( Inf))
@eval typemin{P}(::Type{BigFloat{P}}) = $(BigFloat{P}(-Inf))

function nextfloat{P}(x::BigFloat{P})
    z = BigFloat{P}()
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Int32),
          &z, &x, ROUNDING_MODE[end])
    ccall((:mpfr_nextabove, :libmpfr), Int32, (Ptr{BigFloat{P}},), &z) != 0
    return z
end

function prevfloat{P}(x::BigFloat{P})
    z = BigFloat{P}()
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Int32),
          &z, &x, ROUNDING_MODE[end])
    ccall((:mpfr_nextbelow, :libmpfr), Int32, (Ptr{BigFloat{P}},), &z) != 0
    return z
end

eps(::Type{BigFloat}) = nextfloat(BigFloat(1)) - BigFloat(1)
eps{P}(::Type{BigFloat{P}}) = nextfloat(BigFloat{P}(1)) - BigFloat{P}(1)

realmin(::Type{BigFloat}) = nextfloat(zero(BigFloat))
realmax(::Type{BigFloat}) = prevfloat(BigFloat(Inf))
realmin{P}(::Type{BigFloat}) = nextfloat(zero(BigFloat{P}))
realmax{P}(::Type{BigFloat}) = prevfloat(BigFloat{P}(Inf))


function max{P}(x::BigFloat{P}, y::BigFloat{P})
    z = BigFloat{P}()
    ccall((:mpfr_max, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end

function min{P}(x::BigFloat{P}, y::BigFloat{P})
    z = BigFloat{P}()
    ccall((:mpfr_min, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end
