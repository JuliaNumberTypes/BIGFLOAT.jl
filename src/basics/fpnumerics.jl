
widen(::Type{Float64}) = BigFP
widen(::Type{BigFP}) = BigFP
widen{P}(::Type{BigFP{P}}) = BigFP{2*P+24}


maxintfloat(::Type{BigFP}) = BigFP(2)^precision(BigFP)
maxintfloat{P}(x::BigFP{P}) = BigFP{P}(2)^precision(x)

signbit{P}(x::BigFP{P}) = ccall((:mpfr_signbit, :libmpfr), Int32, (Ptr{BigFP{P}},), &x) != 0

function copysign{P}(x::BigFP{P}, y::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_copysign, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end

function exponent{P}(x::BigFP{P})
    if x == 0 || !isfinite(x)
        throw(DomainError())
    end
    # The '- 1' is to make it work as Base.exponent
    return ccall((:mpfr_get_exp, :libmpfr), Clong, (Ptr{BigFP{P}},), &x) - 1
end

function frexp{P}(x::BigFP{P})
    z = BigFP{P}()
    c = Clong[0]
    ccall((:mpfr_frexp, :libmpfr), Int32, (Ptr{Clong}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Cint), c, &z, &x, ROUNDING_MODE[end])
    return (z, c[1])
end

function significand{P}(x::BigFP{P})
    z = BigFP{P}()
    c = Clong[0]
    ccall((:mpfr_frexp, :libmpfr), Int32, (Ptr{Clong}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Cint), c, &z, &x, ROUNDING_MODE[end])
    # Double the significand to make it work as Base.significand
    ccall((:mpfr_mul_si, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Clong, Int32), &z, &z, 2, ROUNDING_MODE[end])
    return z
end


function isinteger{P}(x::BigFP{P})
    return ccall((:mpfr_integer_p, :libmpfr), Int32, (Ptr{BigFP{P}},), &x) != 0
end


function isinf{P}(x::BigFP{P})
    return ccall((:mpfr_inf_p, :libmpfr), Int32, (Ptr{BigFP{P}},), &x) != 0
end

function isnan{P}(x::BigFP{P})
    return ccall((:mpfr_nan_p, :libmpfr), Int32, (Ptr{BigFP{P}},), &x) != 0
end

isfinite{P}(x::BigFP{P}) = !isinf(x) && !isnan(x)

@eval typemax(::Type{BigFP}) = $(BigFP( Inf))
@eval typemin(::Type{BigFP}) = $(BigFP(-Inf))
@eval typemax{P}(::Type{BigFP{P}}) = $(BigFP{P}( Inf))
@eval typemin{P}(::Type{BigFP{P}}) = $(BigFP{P}(-Inf))

function nextfloat{P}(x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32),
          &z, &x, ROUNDING_MODE[end])
    ccall((:mpfr_nextabove, :libmpfr), Int32, (Ptr{BigFP{P}},), &z) != 0
    return z
end

function prevfloat{P}(x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32),
          &z, &x, ROUNDING_MODE[end])
    ccall((:mpfr_nextbelow, :libmpfr), Int32, (Ptr{BigFP{P}},), &z) != 0
    return z
end

eps(::Type{BigFP}) = nextfloat(BigFP(1)) - BigFP(1)
eps{P}(::Type{BigFP{P}}) = nextfloat(BigFP{P}(1)) - BigFP{P}(1)

realmin(::Type{BigFP}) = nextfloat(zero(BigFP))
realmax(::Type{BigFP}) = prevfloat(BigFP(Inf))
realmin{P}(::Type{BigFP}) = nextfloat(zero(BigFP{P}))
realmax{P}(::Type{BigFP}) = prevfloat(BigFP{P}(Inf))


function max{P}(x::BigFP{P}, y::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_max, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end

function min{P}(x::BigFP{P}, y::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_min, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end
