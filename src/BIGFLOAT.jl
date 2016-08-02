# This file is a part of Julia. License is MIT: http://julialang.org/license

module BIGFLOAT

export
    BigFP,
    setprecision,
    big_str

import
    Base: (*), +, -, /, <, <=, ==, >, >=, ^, besselj, besselj0, besselj1, bessely,
        bessely0, bessely1, ceil, cmp, convert, copysign, div,
        exp, exp2, exponent, factorial, floor, fma, hypot, isinteger,
        isfinite, isinf, isnan, ldexp, log, log2, log10, max, min, mod, modf,
        nextfloat, prevfloat, promote_rule, rem, round, show,
        sum, sqrt, string, print, trunc, precision, exp10, expm1,
        gamma, lgamma, digamma, erf, erfc, zeta, eta, log1p, airyai,
        eps, signbit, sin, cos, tan, sec, csc, cot, acos, asin, atan,
        cosh, sinh, tanh, sech, csch, coth, acosh, asinh, atanh, atan2,
        cbrt, typemax, typemin, unsafe_trunc, realmin, realmax, rounding,
        setrounding, maxintfloat, widen, significand, frexp, tryparse

import Base.Rounding: rounding_raw, setrounding_raw

import Base.GMP: ClongMax, CulongMax, CdoubleMax, Limb

import Base.Math.lgamma_r

function __init__()
    try
        # set exponent to full range by default
        set_emin!(get_emin_min())
        set_emax!(get_emax_max())
    catch ex
        Base.showerror_nostdio(ex,
            "WARNING: Error during initialization of module MPFR")
    end
end

const ROUNDING_MODE = Ref{Cint}(0)
const DEFAULT_PRECISION = [256]

# Basic type and initialization definitions

type BigFP{P} <: AbstractFloat
    prec::Clong
    sign::Cint
    exp::Clong
    d::Ptr{Limb}
    function BigFP()
        N = precision(BigFP)
        z = new(zero(Clong), zero(Cint), zero(Clong), C_NULL)
        ccall((:mpfr_init2,:libmpfr), Void, (Ptr{BigFP}, Clong), &z, N)
        finalizer(z, cglobal((:mpfr_clear, :libmpfr)))
        return z
    end
    # Not recommended for general use
    function BigFP(prec::Clong, sign::Cint, exp::Clong, d::Ptr{Void})
        new(prec, sign, exp, d)
    end
end

BigFP{P}(::Type{BigFP{P}}) = BigFP(((Clong)P), zero(Cint), zero(Clong), C_NULL)

precision{P}(::Type{BigFP{P}}) = P
precision{P}(x::BigFP{P}) = P
precision{T<:BigFP}(::Type{T}) = DEFAULT_PRECISION[1]
precision{T<:BigFP}(x::T) = precision(T)



widen(::Type{Float64}) = BigFP{precision{BigFP}}
widen(::Type{BigFP}) = BigFP
widen{P}(::Type{BigFP{P}}) = BigFP{P}

convert{P}(::Type{BigFP{P}}, x::BigFP{P}) = x
convert(::Type{BigFP}, x::BigFP) = x

# convert to BigFP
for (fJ, fC) in ((:si,:Clong), (:ui,:Culong), (:d,:Float64))
    @eval begin
        function convert{P}(::Type{BigFP}, x::($fC))
            z = BigFP{P}(BigFP{P})
            ccall(($(string(:mpfr_set_,fJ)), :libmpfr), Int32, (Ptr{BigFP{P}}, ($fC), Int32), &z, x, ROUNDING_MODE[])
            return z
        end
    end
end

function convert(::Type{BigFP}, x::BigInt)
    z = BigFP()
    ccall((:mpfr_set_z, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigInt}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

convert(::Type{BigFP}, x::Integer) = BigFP(BigInt(x))

convert(::Type{BigFP}, x::Union{Bool,Int8,Int16,Int32}) = BigFP(convert(Clong,x))
convert(::Type{BigFP}, x::Union{UInt8,UInt16,UInt32}) = BigFP(convert(Culong,x))

convert(::Type{BigFP}, x::Union{Float16,Float32}) = BigFP(Float64(x))
convert(::Type{BigFP}, x::Rational) = BigFP(num(x)) / BigFP(den(x))

function tryparse(::Type{BigFP}, s::AbstractString, base::Int=0)
    z = BigFP()
    err = ccall((:mpfr_set_str, :libmpfr), Int32, (Ptr{BigFP}, Cstring, Int32, Int32), &z, s, base, ROUNDING_MODE[])
    err == 0 ? Nullable(z) : Nullable{BigFP}()
end

convert(::Type{Rational}, x::BigFP) = convert(Rational{BigInt}, x)
convert(::Type{AbstractFloat}, x::BigInt) = BigFP(x)

## BigFP -> Integer
function unsafe_cast(::Type{Int64}, x::BigFP, ri::Cint)
    ccall((:__gmpfr_mpfr_get_sj,:libmpfr), Cintmax_t,
          (Ptr{BigFP}, Cint), &x, ri)
end
function unsafe_cast(::Type{UInt64}, x::BigFP, ri::Cint)
    ccall((:__gmpfr_mpfr_get_uj,:libmpfr), Cuintmax_t,
          (Ptr{BigFP}, Cint), &x, ri)
end

function unsafe_cast{T<:Signed}(::Type{T}, x::BigFP, ri::Cint)
    unsafe_cast(Int64, x, ri) % T
end
function unsafe_cast{T<:Unsigned}(::Type{T}, x::BigFP, ri::Cint)
    unsafe_cast(UInt64, x, ri) % T
end

function unsafe_cast(::Type{BigInt}, x::BigFP, ri::Cint)
    # actually safe, just keep naming consistent
    z = BigInt()
    ccall((:mpfr_get_z, :libmpfr), Int32, (Ptr{BigInt}, Ptr{BigFP}, Int32),
          &z, &x, ri)
    z
end
unsafe_cast(::Type{Int128}, x::BigFP, ri::Cint) = Int128(unsafe_cast(BigInt,x,ri))
unsafe_cast(::Type{UInt128}, x::BigFP, ri::Cint) = UInt128(unsafe_cast(BigInt,x,ri))
unsafe_cast{T<:Integer}(::Type{T}, x::BigFP, r::RoundingMode) = unsafe_cast(T,x,to_mpfr(r))

unsafe_trunc{T<:Integer}(::Type{T}, x::BigFP) = unsafe_cast(T,x,RoundToZero)

function trunc{T<:Union{Signed,Unsigned}}(::Type{T}, x::BigFP)
    (typemin(T) <= x <= typemax(T)) || throw(InexactError())
    unsafe_cast(T,x,RoundToZero)
end
function floor{T<:Union{Signed,Unsigned}}(::Type{T}, x::BigFP)
    (typemin(T) <= x <= typemax(T)) || throw(InexactError())
    unsafe_cast(T,x,RoundDown)
end
function ceil{T<:Union{Signed,Unsigned}}(::Type{T}, x::BigFP)
    (typemin(T) <= x <= typemax(T)) || throw(InexactError())
    unsafe_cast(T,x,RoundUp)
end

function round{T<:Union{Signed,Unsigned}}(::Type{T}, x::BigFP)
    (typemin(T) <= x <= typemax(T)) || throw(InexactError())
    unsafe_cast(T,x,ROUNDING_MODE[])
end

trunc(::Type{BigInt}, x::BigFP) = unsafe_cast(BigInt, x, RoundToZero)
floor(::Type{BigInt}, x::BigFP) = unsafe_cast(BigInt, x, RoundDown)
ceil(::Type{BigInt}, x::BigFP) = unsafe_cast(BigInt, x, RoundUp)
round(::Type{BigInt}, x::BigFP) = unsafe_cast(BigInt, x, ROUNDING_MODE[])

# convert/round/trunc/floor/ceil(Integer, x) should return a BigInt
trunc(::Type{Integer}, x::BigFP) = trunc(BigInt, x)
floor(::Type{Integer}, x::BigFP) = floor(BigInt, x)
ceil(::Type{Integer}, x::BigFP) = ceil(BigInt, x)
round(::Type{Integer}, x::BigFP) = round(BigInt, x)

convert(::Type{Bool}, x::BigFP) = x==0 ? false : x==1 ? true : throw(InexactError())
function convert(::Type{BigInt},x::BigFP)
    isinteger(x) || throw(InexactError())
    trunc(BigInt,x)
end

function convert{T<:Integer}(::Type{T},x::BigFP)
    isinteger(x) || throw(InexactError())
    trunc(T,x)
end

## BigFP -> AbstractFloat
convert(::Type{Float64}, x::BigFP) =
    ccall((:mpfr_get_d,:libmpfr), Float64, (Ptr{BigFP},Int32), &x, ROUNDING_MODE[])
convert(::Type{Float32}, x::BigFP) =
    ccall((:mpfr_get_flt,:libmpfr), Float32, (Ptr{BigFP},Int32), &x, ROUNDING_MODE[])
# TODO: avoid double rounding
convert(::Type{Float16}, x::BigFP) = convert(Float16, convert(Float32, x))

(::Type{Float64})(x::BigFP, r::RoundingMode) =
    ccall((:mpfr_get_d,:libmpfr), Float64, (Ptr{BigFP},Int32), &x, to_mpfr(r))
(::Type{Float32})(x::BigFP, r::RoundingMode) =
    ccall((:mpfr_get_flt,:libmpfr), Float32, (Ptr{BigFP},Int32), &x, to_mpfr(r))
# TODO: avoid double rounding
(::Type{Float16})(x::BigFP, r::RoundingMode) =
    convert(Float16, Float32(x, r))

promote_rule{T<:Real}(::Type{BigFP}, ::Type{T}) = BigFP
promote_rule{T<:AbstractFloat}(::Type{BigInt},::Type{T}) = BigFP
promote_rule{T<:AbstractFloat}(::Type{BigFP},::Type{T}) = BigFP

function convert(::Type{Rational{BigInt}}, x::AbstractFloat)
    if isnan(x); return zero(BigInt)//zero(BigInt); end
    if isinf(x); return copysign(one(BigInt),x)//zero(BigInt); end
    if x == 0;   return zero(BigInt) // one(BigInt); end
    s = max(precision(x) - exponent(x), 0)
    BigInt(ldexp(x,s)) // (BigInt(1) << s)
end

# Basic arithmetic without promotion
for (fJ, fC) in ((:+,:add), (:*,:mul))
    @eval begin
        # BigFP
        function ($fJ)(x::BigFP, y::BigFP)
            z = BigFP()
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, &y, ROUNDING_MODE[])
            return z
        end

        # Unsigned Integer
        function ($fJ)(x::BigFP, c::CulongMax)
            z = BigFP()
            ccall(($(string(:mpfr_,fC,:_ui)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Culong, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        ($fJ)(c::CulongMax, x::BigFP) = ($fJ)(x,c)

        # Signed Integer
        function ($fJ)(x::BigFP, c::ClongMax)
            z = BigFP()
            ccall(($(string(:mpfr_,fC,:_si)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Clong, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        ($fJ)(c::ClongMax, x::BigFP) = ($fJ)(x,c)

        # Float32/Float64
        function ($fJ)(x::BigFP, c::CdoubleMax)
            z = BigFP()
            ccall(($(string(:mpfr_,fC,:_d)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Cdouble, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        ($fJ)(c::CdoubleMax, x::BigFP) = ($fJ)(x,c)

        # BigInt
        function ($fJ)(x::BigFP, c::BigInt)
            z = BigFP()
            ccall(($(string(:mpfr_,fC,:_z)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigInt}, Int32), &z, &x, &c, ROUNDING_MODE[])
            return z
        end
        ($fJ)(c::BigInt, x::BigFP) = ($fJ)(x,c)
    end
end

for (fJ, fC) in ((:-,:sub), (:/,:div))
    @eval begin
        # BigFP
        function ($fJ)(x::BigFP, y::BigFP)
            z = BigFP()
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, &y, ROUNDING_MODE[])
            return z
        end

        # Unsigned Int
        function ($fJ)(x::BigFP, c::CulongMax)
            z = BigFP()
            ccall(($(string(:mpfr_,fC,:_ui)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Culong, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        function ($fJ)(c::CulongMax, x::BigFP)
            z = BigFP()
            ccall(($(string(:mpfr_,:ui_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Culong, Ptr{BigFP}, Int32), &z, c, &x, ROUNDING_MODE[])
            return z
        end

        # Signed Integer
        function ($fJ)(x::BigFP, c::ClongMax)
            z = BigFP()
            ccall(($(string(:mpfr_,fC,:_si)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Clong, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        function ($fJ)(c::ClongMax, x::BigFP)
            z = BigFP()
            ccall(($(string(:mpfr_,:si_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Clong, Ptr{BigFP}, Int32), &z, c, &x, ROUNDING_MODE[])
            return z
        end

        # Float32/Float64
        function ($fJ)(x::BigFP, c::CdoubleMax)
            z = BigFP()
            ccall(($(string(:mpfr_,fC,:_d)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Cdouble, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        function ($fJ)(c::CdoubleMax, x::BigFP)
            z = BigFP()
            ccall(($(string(:mpfr_,:d_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Cdouble, Ptr{BigFP}, Int32), &z, c, &x, ROUNDING_MODE[])
            return z
        end

        # BigInt
        function ($fJ)(x::BigFP, c::BigInt)
            z = BigFP()
            ccall(($(string(:mpfr_,fC,:_z)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigInt}, Int32), &z, &x, &c, ROUNDING_MODE[])
            return z
        end
        # no :mpfr_z_div function
    end
end

function -(c::BigInt, x::BigFP)
    z = BigFP()
    ccall((:mpfr_z_sub, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigInt}, Ptr{BigFP}, Int32), &z, &c, &x, ROUNDING_MODE[])
    return z
end

function fma(x::BigFP, y::BigFP, z::BigFP)
    r = BigFP()
    ccall(("mpfr_fma",:libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &r, &x, &y, &z, ROUNDING_MODE[])
    return r
end

# div
# BigFP
function div(x::BigFP, y::BigFP)
    z = BigFP()
    ccall((:mpfr_div,:libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, &y, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &z)
    return z
end

# Unsigned Int
function div(x::BigFP, c::CulongMax)
    z = BigFP()
    ccall((:mpfr_div_ui, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Culong, Int32), &z, &x, c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &z)
    return z
end
function div(c::CulongMax, x::BigFP)
    z = BigFP()
    ccall((:mpfr_ui_div, :libmpfr), Int32, (Ptr{BigFP}, Culong, Ptr{BigFP}, Int32), &z, c, &x, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &z)
    return z
end

# Signed Integer
function div(x::BigFP, c::ClongMax)
    z = BigFP()
    ccall((:mpfr_div_si, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Clong, Int32), &z, &x, c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &z)
    return z
end
function div(c::ClongMax, x::BigFP)
    z = BigFP()
    ccall((:mpfr_si_div, :libmpfr), Int32, (Ptr{BigFP}, Clong, Ptr{BigFP}, Int32), &z, c, &x, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &z)
    return z
end

# Float32/Float64
function div(x::BigFP, c::CdoubleMax)
    z = BigFP()
    ccall((:mpfr_div_d, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Cdouble, Int32), &z, &x, c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &z)
    return z
end
function div(c::CdoubleMax, x::BigFP)
    z = BigFP()
    ccall((:mpfr_d_div, :libmpfr), Int32, (Ptr{BigFP}, Cdouble, Ptr{BigFP}, Int32), &z, c, &x, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &z)
    return z
end

# BigInt
function div(x::BigFP, c::BigInt)
    z = BigFP()
    ccall((:mpfr_div_z, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigInt}, Int32), &z, &x, &c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &z)
    return z
end


# More efficient commutative operations
for (fJ, fC, fI) in ((:+, :add, 0), (:*, :mul, 1))
    @eval begin
        function ($fJ)(a::BigFP, b::BigFP, c::BigFP)
            z = BigFP()
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &a, &b, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &z, &c, ROUNDING_MODE[])
            return z
        end
        function ($fJ)(a::BigFP, b::BigFP, c::BigFP, d::BigFP)
            z = BigFP()
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &a, &b, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &z, &c, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &z, &d, ROUNDING_MODE[])
            return z
        end
        function ($fJ)(a::BigFP, b::BigFP, c::BigFP, d::BigFP, e::BigFP)
            z = BigFP()
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &a, &b, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &z, &c, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &z, &d, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &z, &e, ROUNDING_MODE[])
            return z
        end
    end
end

function -(x::BigFP)
    z = BigFP()
    ccall((:mpfr_neg, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function sqrt(x::BigFP)
    isnan(x) && return x
    z = BigFP()
    ccall((:mpfr_sqrt, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
    if isnan(z)
        throw(DomainError())
    end
    return z
end

sqrt(x::BigInt) = sqrt(BigFP(x))

function ^(x::BigFP, y::BigFP)
    z = BigFP()
    ccall((:mpfr_pow, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function ^(x::BigFP, y::CulongMax)
    z = BigFP()
    ccall((:mpfr_pow_ui, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Culong, Int32), &z, &x, y, ROUNDING_MODE[])
    return z
end

function ^(x::BigFP, y::ClongMax)
    z = BigFP()
    ccall((:mpfr_pow_si, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Clong, Int32), &z, &x, y, ROUNDING_MODE[])
    return z
end

function ^(x::BigFP, y::BigInt)
    z = BigFP()
    ccall((:mpfr_pow_z, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigInt}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

^(x::BigFP, y::Integer)  = typemin(Clong)  <= y <= typemax(Clong)  ? x^Clong(y)  : x^BigInt(y)
^(x::BigFP, y::Unsigned) = typemin(Culong) <= y <= typemax(Culong) ? x^Culong(y) : x^BigInt(y)

for f in (:exp, :exp2, :exp10, :expm1, :digamma, :erf, :erfc, :zeta,
          :cosh,:sinh,:tanh,:sech,:csch,:coth, :cbrt)
    @eval function $f(x::BigFP)
        z = BigFP()
        ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
        return z
    end
end

# return log(2)
function big_ln2()
    c = BigFP()
    ccall((:mpfr_const_log2, :libmpfr), Cint, (Ptr{BigFP}, Int32),
          &c, MPFR.ROUNDING_MODE[])
    return c
end

function eta(x::BigFP)
    x == 1 && return big_ln2()
    return -zeta(x) * expm1(big_ln2()*(1-x))
end

function airyai(x::BigFP)
    z = BigFP()
    ccall((:mpfr_ai, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end
airy(x::BigFP) = airyai(x)

function ldexp(x::BigFP, n::Clong)
    z = BigFP()
    ccall((:mpfr_mul_2si, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Clong, Int32), &z, &x, n, ROUNDING_MODE[])
    return z
end
function ldexp(x::BigFP, n::Culong)
    z = BigFP()
    ccall((:mpfr_mul_2ui, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Culong, Int32), &z, &x, n, ROUNDING_MODE[])
    return z
end
ldexp(x::BigFP, n::ClongMax) = ldexp(x, convert(Clong, n))
ldexp(x::BigFP, n::CulongMax) = ldexp(x, convert(Culong, n))
ldexp(x::BigFP, n::Integer) = x*exp2(BigFP(n))

function besselj0(x::BigFP)
    z = BigFP()
    ccall((:mpfr_j0, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function besselj1(x::BigFP)
    z = BigFP()
    ccall((:mpfr_j1, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function besselj(n::Integer, x::BigFP)
    z = BigFP()
    ccall((:mpfr_jn, :libmpfr), Int32, (Ptr{BigFP}, Clong, Ptr{BigFP}, Int32), &z, n, &x, ROUNDING_MODE[])
    return z
end

function bessely0(x::BigFP)
    if x < 0
        throw(DomainError())
    end
    z = BigFP()
    ccall((:mpfr_y0, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function bessely1(x::BigFP)
    if x < 0
        throw(DomainError())
    end
    z = BigFP()
    ccall((:mpfr_y1, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function bessely(n::Integer, x::BigFP)
    if x < 0
        throw(DomainError())
    end
    z = BigFP()
    ccall((:mpfr_yn, :libmpfr), Int32, (Ptr{BigFP}, Clong, Ptr{BigFP}, Int32), &z, n, &x, ROUNDING_MODE[])
    return z
end

function factorial(x::BigFP)
    if x < 0 || !isinteger(x)
        throw(DomainError())
    end
    ui = convert(Culong, x)
    z = BigFP()
    ccall((:mpfr_fac_ui, :libmpfr), Int32, (Ptr{BigFP}, Culong, Int32), &z, ui, ROUNDING_MODE[])
    return z
end

function hypot(x::BigFP, y::BigFP)
    z = BigFP()
    ccall((:mpfr_hypot, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

for f in (:log, :log2, :log10)
    @eval function $f(x::BigFP)
        if x < 0
            throw(DomainError())
        end
        z = BigFP()
        ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
        return z
    end
end

function log1p(x::BigFP)
    if x < -1
        throw(DomainError())
    end
    z = BigFP()
    ccall((:mpfr_log1p, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function max(x::BigFP, y::BigFP)
    z = BigFP()
    ccall((:mpfr_max, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function min(x::BigFP, y::BigFP)
    z = BigFP()
    ccall((:mpfr_min, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function modf(x::BigFP)
    if isinf(x)
        return (BigFP(NaN), x)
    end
    zint = BigFP()
    zfloat = BigFP()
    ccall((:mpfr_modf, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &zint, &zfloat, &x, ROUNDING_MODE[])
    return (zfloat, zint)
end

function rem(x::BigFP, y::BigFP)
    z = BigFP()
    ccall((:mpfr_fmod, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function sum(arr::AbstractArray{BigFP})
    z = BigFP(0)
    for i in arr
        ccall((:mpfr_add, :libmpfr), Int32,
            (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Cint),
            &z, &z, &i, 0)
    end
    return z
end

# Functions for which NaN results are converted to DomainError, following Base
for f in (:sin,:cos,:tan,:sec,:csc,
          :acos,:asin,:atan,:acosh,:asinh,:atanh, :gamma)
    @eval begin
        function ($f)(x::BigFP)
            if isnan(x)
                return x
            end
            z = BigFP()
            ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, ROUNDING_MODE[])
            if isnan(z)
                throw(DomainError())
            end
            return z
        end
    end
end

# log of absolute value of gamma function
const lgamma_signp = Array{Cint}(1)
function lgamma(x::BigFP)
    z = BigFP()
    ccall((:mpfr_lgamma,:libmpfr), Cint, (Ptr{BigFP}, Ptr{Cint}, Ptr{BigFP}, Int32), &z, lgamma_signp, &x, ROUNDING_MODE[])
    return z
end

lgamma_r(x::BigFP) = (lgamma(x), lgamma_signp[1])

function atan2(y::BigFP, x::BigFP)
    z = BigFP()
    ccall((:mpfr_atan2, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &y, &x, ROUNDING_MODE[])
    return z
end

# Utility functions
==(x::BigFP, y::BigFP) = ccall((:mpfr_equal_p, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &x, &y) != 0
<=(x::BigFP, y::BigFP) = ccall((:mpfr_lessequal_p, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &x, &y) != 0
>=(x::BigFP, y::BigFP) = ccall((:mpfr_greaterequal_p, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &x, &y) != 0
<(x::BigFP, y::BigFP) = ccall((:mpfr_less_p, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &x, &y) != 0
>(x::BigFP, y::BigFP) = ccall((:mpfr_greater_p, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &x, &y) != 0

function cmp(x::BigFP, y::BigInt)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_z, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigInt}), &x, &y)
end
function cmp(x::BigFP, y::ClongMax)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_si, :libmpfr), Int32, (Ptr{BigFP}, Clong), &x, y)
end
function cmp(x::BigFP, y::CulongMax)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_ui, :libmpfr), Int32, (Ptr{BigFP}, Culong), &x, y)
end
cmp(x::BigFP, y::Integer) = cmp(x,big(y))
cmp(x::Integer, y::BigFP) = -cmp(y,x)

function cmp(x::BigFP, y::CdoubleMax)
    (isnan(x) || isnan(y)) && throw(DomainError())
    ccall((:mpfr_cmp_d, :libmpfr), Int32, (Ptr{BigFP}, Cdouble), &x, y)
end
cmp(x::CdoubleMax, y::BigFP) = -cmp(y,x)

==(x::BigFP, y::Integer)   = !isnan(x) && cmp(x,y) == 0
==(x::Integer, y::BigFP)   = y == x
==(x::BigFP, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) == 0
==(x::CdoubleMax, y::BigFP) = y == x

<(x::BigFP, y::Integer)   = !isnan(x) && cmp(x,y) < 0
<(x::Integer, y::BigFP)   = !isnan(y) && cmp(y,x) > 0
<(x::BigFP, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) < 0
<(x::CdoubleMax, y::BigFP) = !isnan(x) && !isnan(y) && cmp(y,x) > 0

<=(x::BigFP, y::Integer)   = !isnan(x) && cmp(x,y) <= 0
<=(x::Integer, y::BigFP)   = !isnan(y) && cmp(y,x) >= 0
<=(x::BigFP, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) <= 0
<=(x::CdoubleMax, y::BigFP) = !isnan(x) && !isnan(y) && cmp(y,x) >= 0

signbit(x::BigFP) = ccall((:mpfr_signbit, :libmpfr), Int32, (Ptr{BigFP},), &x) != 0

function precision(x::BigFP)  # precision of an object of type BigFP
    return ccall((:mpfr_get_prec, :libmpfr), Clong, (Ptr{BigFP},), &x)
end

precision(::Type{BigFP}) = DEFAULT_PRECISION[end]  # precision of the type BigFP itself

"""
    setprecision([T=BigFP,] precision::Int)

Set the precision (in bits) to be used for `T` arithmetic.
"""
function setprecision(::Type{BigFP}, precision::Int)
    if precision < 2
        throw(DomainError())
    end
    DEFAULT_PRECISION[end] = precision
end

setprecision(precision::Int) = setprecision(BigFP, precision)

maxintfloat(x::BigFP) = BigFP(2)^precision(x)
maxintfloat(::Type{BigFP}) = BigFP(2)^precision(BigFP)

to_mpfr(::RoundingMode{:Nearest}) = Cint(0)
to_mpfr(::RoundingMode{:ToZero}) = Cint(1)
to_mpfr(::RoundingMode{:Up}) = Cint(2)
to_mpfr(::RoundingMode{:Down}) = Cint(3)
to_mpfr(::RoundingMode{:FromZero}) = Cint(4)

function from_mpfr(c::Integer)
    if c == 0
        return RoundNearest
    elseif c == 1
        return RoundToZero
    elseif c == 2
        return RoundUp
    elseif c == 3
        return RoundDown
    elseif c == 4
        return RoundFromZero
    else
        throw(ArgumentError("invalid MPFR rounding mode code: $c"))
    end
    RoundingMode(c)
end

rounding_raw(::Type{BigFP}) = ROUNDING_MODE[]
setrounding_raw(::Type{BigFP},i::Integer) = ROUNDING_MODE[] = i

rounding(::Type{BigFP}) = from_mpfr(rounding_raw(BigFP))
setrounding(::Type{BigFP},r::RoundingMode) = setrounding_raw(BigFP,to_mpfr(r))

function copysign(x::BigFP, y::BigFP)
    z = BigFP()
    ccall((:mpfr_copysign, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Ptr{BigFP}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function exponent(x::BigFP)
    if x == 0 || !isfinite(x)
        throw(DomainError())
    end
    # The '- 1' is to make it work as Base.exponent
    return ccall((:mpfr_get_exp, :libmpfr), Clong, (Ptr{BigFP},), &x) - 1
end

function frexp(x::BigFP)
    z = BigFP()
    c = Ref{Clong}()
    ccall((:mpfr_frexp, :libmpfr), Int32, (Ptr{Clong}, Ptr{BigFP}, Ptr{BigFP}, Cint), c, &z, &x, ROUNDING_MODE[])
    return (z, c[])
end

function significand(x::BigFP)
    z = BigFP()
    c = Ref{Clong}()
    ccall((:mpfr_frexp, :libmpfr), Int32, (Ptr{Clong}, Ptr{BigFP}, Ptr{BigFP}, Cint), c, &z, &x, ROUNDING_MODE[])
    # Double the significand to make it work as Base.significand
    ccall((:mpfr_mul_si, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Clong, Int32), &z, &z, 2, ROUNDING_MODE[])
    return z
end

function isinteger(x::BigFP)
    return ccall((:mpfr_integer_p, :libmpfr), Int32, (Ptr{BigFP},), &x) != 0
end

for f in (:ceil, :floor, :trunc)
    @eval begin
        function ($f)(x::BigFP)
            z = BigFP()
            ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &x)
            return z
        end
    end
end

function round(x::BigFP)
    z = BigFP()
    ccall((:mpfr_rint, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Cint), &z, &x, ROUNDING_MODE[])
    return z
end
function round(x::BigFP,::RoundingMode{:NearestTiesAway})
    z = BigFP()
    ccall((:mpfr_round, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}), &z, &x)
    return z
end

function isinf(x::BigFP)
    return ccall((:mpfr_inf_p, :libmpfr), Int32, (Ptr{BigFP},), &x) != 0
end

function isnan(x::BigFP)
    return ccall((:mpfr_nan_p, :libmpfr), Int32, (Ptr{BigFP},), &x) != 0
end

isfinite(x::BigFP) = !isinf(x) && !isnan(x)

@eval typemax(::Type{BigFP}) = $(BigFP( Inf))
@eval typemin(::Type{BigFP}) = $(BigFP(-Inf))

function nextfloat(x::BigFP)
    z = BigFP()
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32),
          &z, &x, ROUNDING_MODE[])
    ccall((:mpfr_nextabove, :libmpfr), Int32, (Ptr{BigFP},), &z) != 0
    return z
end

function prevfloat(x::BigFP)
    z = BigFP()
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32),
          &z, &x, ROUNDING_MODE[])
    ccall((:mpfr_nextbelow, :libmpfr), Int32, (Ptr{BigFP},), &z) != 0
    return z
end

eps(::Type{BigFP}) = nextfloat(BigFP(1)) - BigFP(1)

realmin(::Type{BigFP}) = nextfloat(zero(BigFP))
realmax(::Type{BigFP}) = prevfloat(BigFP(Inf))

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
function setprecision{T}(f::Function, ::Type{T}, prec::Integer)
    old_prec = precision(T)
    setprecision(T, prec)
    try
        return f()
    finally
        setprecision(T, old_prec)
    end
end

setprecision(f::Function, precision::Integer) = setprecision(f, BigFP, precision)

function string(x::BigFP)
    # In general, the number of decimal places needed to read back the number exactly
    # is, excluding the most significant, ceil(log(10, 2^precision(x)))
    k = ceil(Int32, precision(x) * 0.3010299956639812)
    lng = k + Int32(8) # Add space for the sign, the most significand digit, the dot and the exponent
    buf = Array{UInt8}(lng + 1)
    # format strings are guaranteed to contain no NUL, so we don't use Cstring
    lng = ccall((:mpfr_snprintf,:libmpfr), Int32, (Ptr{UInt8}, Culong, Ptr{UInt8}, Ptr{BigFP}...), buf, lng + 1, "%.Re", &x)
    if lng < k + 5 # print at least k decimal places
        lng = ccall((:mpfr_sprintf,:libmpfr), Int32, (Ptr{UInt8}, Ptr{UInt8}, Ptr{BigFP}...), buf, "%.$(k)Re", &x)
    elseif lng > k + 8
        buf = Array{UInt8}(lng + 1)
        lng = ccall((:mpfr_snprintf,:libmpfr), Int32, (Ptr{UInt8}, Culong, Ptr{UInt8}, Ptr{BigFP}...), buf, lng + 1, "%.Re", &x)
    end
    n = (1 <= x < 10 || -10 < x <= -1 || x == 0) ? lng - 4 : lng
    return String(buf[1:n])
end

print(io::IO, b::BigFP) = print(io, string(b))
show(io::IO, b::BigFP) = print(io, string(b))

# get/set exponent min/max
get_emax() = ccall((:mpfr_get_emax, :libmpfr), Clong, ())
get_emax_min() = ccall((:mpfr_get_emax_min, :libmpfr), Clong, ())
get_emax_max() = ccall((:mpfr_get_emax_max, :libmpfr), Clong, ())

get_emin() = ccall((:mpfr_get_emin, :libmpfr), Clong, ())
get_emin_min() = ccall((:mpfr_get_emin_min, :libmpfr), Clong, ())
get_emin_max() = ccall((:mpfr_get_emin_max, :libmpfr), Clong, ())

set_emax!(x) = ccall((:mpfr_set_emax, :libmpfr), Void, (Clong,), x)
set_emin!(x) = ccall((:mpfr_set_emin, :libmpfr), Void, (Clong,), x)

function Base.deepcopy_internal(x::BigFP, stackdict::ObjectIdDict)
    if haskey(stackdict, x)
        return stackdict[x]
    end
    N = precision(x)
    y = BigFP(zero(Clong), zero(Cint), zero(Clong), C_NULL)
    ccall((:mpfr_init2,:libmpfr), Void, (Ptr{BigFP}, Clong), &y, N)
    finalizer(y, cglobal((:mpfr_clear, :libmpfr)))
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32),
          &y, &x, ROUNDING_MODE[])
    stackdict[x] = y
    return y
end

end #module
