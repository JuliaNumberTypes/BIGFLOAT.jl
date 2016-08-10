# This file is a part of Julia. License is MIT: http://julialang.org/license

module MPFR3

export
    BigFloat,
    setprecision,
    big_str

import
    Base: @big_str, *, +, -, /, <, <=, ==, >, >=, ^,
          besselj, besselj0, besselj1, bessely, bessely0, bessely1,
          ceil, cmp, convert, copysign, div, cbrt, typemax, typemin,
          unsafe_trunc, realmin, realmax, rounding,
          cosh, sinh, tanh, sech, csch, coth, acosh, asinh, atanh, atan2,
          eps, signbit, sin, cos, tan, sec, csc, cot, acos, asin, atan,
          exp, exp2, exponent, factorial, floor, fma, hypot, isinteger,
          gamma, lgamma, digamma, erf, erfc, zeta, eta, airyai,
          isfinite, isinf, isnan, ldexp, log, log2, log10, log1p, exp10, expm1,
          max, min, mod, modf, nextfloat, prevfloat, promote_rule, rem, round,
          show, setrounding, maxintfloat, widen, significand, frexp, tryparse,
          sum, sqrt, string, print, trunc, precision, setprecision

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

type BigFloat{P} <: AbstractFloat
    prec::Clong
    sign::Cint
    exp::Clong
    d::Ptr{Limb}
    function BigFloat()
    N = Clong(P)
    z = BigFloat{P}(zero(Clong), zero(Cint), zero(Clong), C_NULL)
    ccall((:mpfr_init2,:libmpfr), Void, (Ptr{BigFloat}, Clong), &z, N)
    finalizer(z, cglobal((:mpfr_clear, :libmpfr)))
    return z
    end
    # Not recommended for general use
    function BigFloat(prec::Clong, sign::Cint, exp::Clong, d::Ptr{Void})
        new(prec, sign, exp, d)
    end
end


precision{P}(::Type{BigFloat{P}}) = P
precision{P}(x::BigFloat{P}) = P

function initializer{P}(::Type{BigFloat{P}})
    N = Clong(P)
    z = BigFloat{P}(zero(Clong), zero(Cint), zero(Clong), C_NULL)
    ccall((:mpfr_init2,:libmpfr), Void, (Ptr{BigFloat}, Clong), &z, N)
    finalizer(z, cglobal((:mpfr_clear, :libmpfr)))
    return z
end
initializer(::Type{BigFloat}) = initializer(BigFloat{precision(BigFloat)})

BigFloat{P}(::Type{BigFloat{P}}) = initializer(BigFloat{P})
BigFloat() = initializer(BigFloat)
function BigFloat(prec::Int64, sign::Int32, exp::Int64, d::Ptr{Void})
    P = Int(prec)
    return BigFloat{P}(prec, sign, exp, d)
end




widen(::Type{Float64}) = BigFloat
widen(::Type{BigFloat}) = BigFloat

convert{BF<:BigFloat}(::Type{BF}, x::BF) = x

# convert to BF
for (fJ, fC) in ((:si,:Clong), (:ui,:Culong), (:d,:Float64))
    @eval begin
        function convert{BF<:BigFloat}(::Type{BF}, x::($fC))
            z = BF()
            ccall(($(string(:mpfr_set_,fJ)), :libmpfr), Int32, (Ptr{BF}, ($fC), Int32), &z, x, ROUNDING_MODE[])
            return z
        end
    end
end

function convert{BF<:BigFloat}(::Type{BF}, x::BigInt)
    z = BF()
    ccall((:mpfr_set_z, :libmpfr), Int32, (Ptr{BF}, Ptr{BigInt}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

convert{BF<:BigFloat}(::Type{BF}, x::Integer) = BF(BigInt(x))

convert{BF<:BigFloat}(::Type{BF}, x::Union{Bool,Int8,Int16,Int32}) = BF(convert(Clong,x))
convert{BF<:BigFloat}(::Type{BF}, x::Union{UInt8,UInt16,UInt32}) = BF(convert(Culong,x))

convert{BF<:BigFloat}(::Type{BF}, x::Union{Float16,Float32}) = BF(Float64(x))
convert{BF<:BigFloat}(::Type{BF}, x::Rational) = BF(num(x)) / BF(den(x))

function tryparse{BF<:BigFloat}(::Type{BF}, s::AbstractString, base::Int=0)
    z = BF()
    err = ccall((:mpfr_set_str, :libmpfr), Int32, (Ptr{BF}, Cstring, Int32, Int32), &z, s, base, ROUNDING_MODE[])
    err == 0 ? Nullable(z) : Nullable{BF}()
end

convert{BF<:BigFloat}(::Type{Rational}, x::BF) = convert(Rational{BigInt}, x)
convert(::Type{AbstractFloat}, x::BigInt) = BigFloat(x)

## BF -> Integer
function unsafe_cast{BF<:BigFloat}(::Type{Int64}, x::BF, ri::Cint)
    ccall((:__gmpfr_mpfr_get_sj,:libmpfr), Cintmax_t,
          (Ptr{BF}, Cint), &x, ri)
end
function unsafe_cast{BF<:BigFloat}(::Type{UInt64}, x::BF, ri::Cint)
    ccall((:__gmpfr_mpfr_get_uj,:libmpfr), Cuintmax_t,
          (Ptr{BF}, Cint), &x, ri)
end

function unsafe_cast{T<:Signed,BF<:BigFloat}(::Type{T}, x::BF, ri::Cint)
    unsafe_cast(Int64, x, ri) % T
end
function unsafe_cast{T<:Unsigned,BF<:BigFloat}(::Type{T}, x::BF, ri::Cint)
    unsafe_cast(UInt64, x, ri) % T
end

function unsafe_cast{BF<:BigFloat}(::Type{BigInt}, x::BF, ri::Cint)
    # actually safe, just keep naming consistent
    z = BigInt()
    ccall((:mpfr_get_z, :libmpfr), Int32, (Ptr{BigInt}, Ptr{BF}, Int32),
          &z, &x, ri)
    z
end
unsafe_cast{BF<:BigFloat}(::Type{Int128}, x::BF, ri::Cint) = Int128(unsafe_cast(BigInt,x,ri))
unsafe_cast{BF<:BigFloat}(::Type{UInt128}, x::BF, ri::Cint) = UInt128(unsafe_cast(BigInt,x,ri))
unsafe_cast{T<:Integer,BF<:BigFloat}(::Type{T}, x::BF, r::RoundingMode) = unsafe_cast(T,x,to_mpfr(r))

unsafe_trunc{T<:Integer,BF<:BigFloat}(::Type{T}, x::BF) = unsafe_cast(T,x,RoundToZero)

function trunc{T<:Union{Signed,Unsigned},BF<:BigFloat}(::Type{T}, x::BF)
    (typemin(T) <= x <= typemax(T)) || throw(InexactError())
    unsafe_cast(T,x,RoundToZero)
end
function floor{T<:Union{Signed,Unsigned},BF<:BigFloat}(::Type{T}, x::BF)
    (typemin(T) <= x <= typemax(T)) || throw(InexactError())
    unsafe_cast(T,x,RoundDown)
end
function ceil{T<:Union{Signed,Unsigned},BF<:BigFloat}(::Type{T}, x::BF)
    (typemin(T) <= x <= typemax(T)) || throw(InexactError())
    unsafe_cast(T,x,RoundUp)
end

function round{T<:Union{Signed,Unsigned},BF<:BigFloat}(::Type{T}, x::BF)
    (typemin(T) <= x <= typemax(T)) || throw(InexactError())
    unsafe_cast(T,x,ROUNDING_MODE[])
end

trunc{BF<:BigFloat}(::Type{BigInt}, x::BF) = unsafe_cast(BigInt, x, RoundToZero)
floor{BF<:BigFloat}(::Type{BigInt}, x::BF) = unsafe_cast(BigInt, x, RoundDown)
ceil{BF<:BigFloat}(::Type{BigInt}, x::BF) = unsafe_cast(BigInt, x, RoundUp)
round{BF<:BigFloat}(::Type{BigInt}, x::BF) = unsafe_cast(BigInt, x, ROUNDING_MODE[])

# convert/round/trunc/floor/ceil(Integer, x) should return a BigInt
trunc{BF<:BigFloat}(::Type{Integer}, x::BF) = trunc(BigInt, x)
floor{BF<:BigFloat}(::Type{Integer}, x::BF) = floor(BigInt, x)
ceil{BF<:BigFloat}(::Type{Integer}, x::BF) = ceil(BigInt, x)
round{BF<:BigFloat}(::Type{Integer}, x::BF) = round(BigInt, x)

convert{BF<:BigFloat}(::Type{Bool}, x::BF) = x==0 ? false : x==1 ? true : throw(InexactError())
function convert{BF<:BigFloat}(::Type{BigInt},x::BF)
    isinteger(x) || throw(InexactError())
    trunc(BigInt,x)
end

function convert{T<:Integer,BF<:BigFloat}(::Type{T},x::BF)
    isinteger(x) || throw(InexactError())
    trunc(T,x)
end

## BF -> AbstractFloat
convert{BF<:BigFloat}(::Type{Float64}, x::BF) =
    ccall((:mpfr_get_d,:libmpfr), Float64, (Ptr{BF},Int32), &x, ROUNDING_MODE[])
convert{BF<:BigFloat}(::Type{Float32}, x::BF) =
    ccall((:mpfr_get_flt,:libmpfr), Float32, (Ptr{BF},Int32), &x, ROUNDING_MODE[])
# TODO: avoid double rounding
convert{BF<:BigFloat}(::Type{Float16}, x::BF) = convert(Float16, convert(Float32, x))

(::Type{Float64}){BF<:BigFloat}(x::BF, r::RoundingMode) =
    ccall((:mpfr_get_d,:libmpfr), Float64, (Ptr{BF},Int32), &x, to_mpfr(r))
(::Type{Float32}){BF<:BigFloat}(x::BF, r::RoundingMode) =
    ccall((:mpfr_get_flt,:libmpfr), Float32, (Ptr{BF},Int32), &x, to_mpfr(r))
# TODO: avoid double rounding
(::Type{Float16}){BF<:BigFloat}(x::BF, r::RoundingMode) =
    convert(Float16, Float32(x, r))



promote_rule(::Type{BigInt},::Type{BigFloat}) = BigFloat
promote_rule{T<:Real}(::Type{BigFloat}, ::Type{T}) = BigFloat
promote_rule{T<:AbstractFloat}(::Type{BigInt},::Type{T}) = BigFloat
promote_rule{T<:AbstractFloat}(::Type{BigFloat},::Type{T}) = BigFloat

promote_rule{P}(::Type{BigInt},::Type{BigFloat{P}}) = BigFloat{P}
promote_rule{BF<:BigFloat}(::Type{BigInt},::Type{BF}) = BF
promote_rule{T<:Real,BF<:BigFloat}(::Type{BF}, ::Type{T}) = BF
promote_rule{T<:AbstractFloat,BF<:BigFloat}(::Type{BF},::Type{T}) = BF

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
        # BF
        function ($fJ){BF<:BigFloat}(x::BF, y::BF)
            z = BF()
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &x, &y, ROUNDING_MODE[])
            return z
        end

        # Unsigned Integer
        function ($fJ){BF<:BigFloat}(x::BF, c::CulongMax)
            z = BF()
            ccall(($(string(:mpfr_,fC,:_ui)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Culong, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        ($fJ){BF<:BigFloat}(c::CulongMax, x::BF) = ($fJ)(x,c)

        # Signed Integer
        function ($fJ){BF<:BigFloat}(x::BF, c::ClongMax)
            z = BF()
            ccall(($(string(:mpfr_,fC,:_si)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Clong, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        ($fJ){BF<:BigFloat}(c::ClongMax, x::BF) = ($fJ)(x,c)

        # Float32/Float64
        function ($fJ){BF<:BigFloat}(x::BF, c::CdoubleMax)
            z = BF()
            ccall(($(string(:mpfr_,fC,:_d)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Cdouble, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        ($fJ){BF<:BigFloat}(c::CdoubleMax, x::BF) = ($fJ)(x,c)

        # BigInt
        function ($fJ){BF<:BigFloat}(x::BF, c::BigInt)
            z = BF()
            ccall(($(string(:mpfr_,fC,:_z)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BigInt}, Int32), &z, &x, &c, ROUNDING_MODE[])
            return z
        end
        ($fJ){BF<:BigFloat}(c::BigInt, x::BF) = ($fJ)(x,c)
    end
end

for (fJ, fC) in ((:-,:sub), (:/,:div))
    @eval begin
        # BF
        function ($fJ){BF<:BigFloat}(x::BF, y::BF)
            z = BF()
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &x, &y, ROUNDING_MODE[])
            return z
        end

        # Unsigned Int
        function ($fJ){BF<:BigFloat}(x::BF, c::CulongMax)
            z = BF()
            ccall(($(string(:mpfr_,fC,:_ui)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Culong, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        function ($fJ){BF<:BigFloat}(c::CulongMax, x::BF)
            z = BF()
            ccall(($(string(:mpfr_,:ui_,fC)), :libmpfr), Int32, (Ptr{BF}, Culong, Ptr{BF}, Int32), &z, c, &x, ROUNDING_MODE[])
            return z
        end

        # Signed Integer
        function ($fJ){BF<:BigFloat}(x::BF, c::ClongMax)
            z = BF()
            ccall(($(string(:mpfr_,fC,:_si)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Clong, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        function ($fJ){BF<:BigFloat}(c::ClongMax, x::BF)
            z = BF()
            ccall(($(string(:mpfr_,:si_,fC)), :libmpfr), Int32, (Ptr{BF}, Clong, Ptr{BF}, Int32), &z, c, &x, ROUNDING_MODE[])
            return z
        end

        # Float32/Float64
        function ($fJ){BF<:BigFloat}(x::BF, c::CdoubleMax)
            z = BF()
            ccall(($(string(:mpfr_,fC,:_d)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Cdouble, Int32), &z, &x, c, ROUNDING_MODE[])
            return z
        end
        function ($fJ){BF<:BigFloat}(c::CdoubleMax, x::BF)
            z = BF()
            ccall(($(string(:mpfr_,:d_,fC)), :libmpfr), Int32, (Ptr{BF}, Cdouble, Ptr{BF}, Int32), &z, c, &x, ROUNDING_MODE[])
            return z
        end

        # BigInt
        function ($fJ){BF<:BigFloat}(x::BF, c::BigInt)
            z = BF()
            ccall(($(string(:mpfr_,fC,:_z)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BigInt}, Int32), &z, &x, &c, ROUNDING_MODE[])
            return z
        end
        # no :mpfr_z_div function
    end
end

function -{BF<:BigFloat}(c::BigInt, x::BF)
    z = BF()
    ccall((:mpfr_z_sub, :libmpfr), Int32, (Ptr{BF}, Ptr{BigInt}, Ptr{BF}, Int32), &z, &c, &x, ROUNDING_MODE[])
    return z
end

function fma{BF<:BigFloat}(x::BF, y::BF, z::BF)
    r = BF()
    ccall(("mpfr_fma",:libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &r, &x, &y, &z, ROUNDING_MODE[])
    return r
end

# div
# BF
function div{BF<:BigFloat}(x::BF, y::BF)
    z = BF()
    ccall((:mpfr_div,:libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &x, &y, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &z)
    return z
end

# Unsigned Int
function div{BF<:BigFloat}(x::BF, c::CulongMax)
    z = BF()
    ccall((:mpfr_div_ui, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Culong, Int32), &z, &x, c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &z)
    return z
end
function div{BF<:BigFloat}(c::CulongMax, x::BF)
    z = BF()
    ccall((:mpfr_ui_div, :libmpfr), Int32, (Ptr{BF}, Culong, Ptr{BF}, Int32), &z, c, &x, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &z)
    return z
end

# Signed Integer
function div{BF<:BigFloat}(x::BF, c::ClongMax)
    z = BF()
    ccall((:mpfr_div_si, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Clong, Int32), &z, &x, c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &z)
    return z
end
function div{BF<:BigFloat}(c::ClongMax, x::BF)
    z = BF()
    ccall((:mpfr_si_div, :libmpfr), Int32, (Ptr{BF}, Clong, Ptr{BF}, Int32), &z, c, &x, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &z)
    return z
end

# Float32/Float64
function div{BF<:BigFloat}(x::BF, c::CdoubleMax)
    z = BF()
    ccall((:mpfr_div_d, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Cdouble, Int32), &z, &x, c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &z)
    return z
end
function div{BF<:BigFloat}(c::CdoubleMax, x::BF)
    z = BF()
    ccall((:mpfr_d_div, :libmpfr), Int32, (Ptr{BF}, Cdouble, Ptr{BF}, Int32), &z, c, &x, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &z)
    return z
end

# BigInt
function div{BF<:BigFloat}(x::BF, c::BigInt)
    z = BF()
    ccall((:mpfr_div_z, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BigInt}, Int32), &z, &x, &c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &z)
    return z
end


# More efficient commutative operations
for (fJ, fC, fI) in ((:+, :add, 0), (:*, :mul, 1))
    @eval begin
        function ($fJ){BF<:BigFloat}(a::BF, b::BF, c::BF)
            z = BF()
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &a, &b, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &z, &c, ROUNDING_MODE[])
            return z
        end
        function ($fJ){BF<:BigFloat}(a::BF, b::BF, c::BF, d::BF)
            z = BF()
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &a, &b, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &z, &c, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &z, &d, ROUNDING_MODE[])
            return z
        end
        function ($fJ){BF<:BigFloat}(a::BF, b::BF, c::BF, d::BF, e::BF)
            z = BF()
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &a, &b, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &z, &c, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &z, &d, ROUNDING_MODE[])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &z, &e, ROUNDING_MODE[])
            return z
        end
    end
end

function -{BF<:BigFloat}(x::BF)
    z = BF()
    ccall((:mpfr_neg, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function sqrt{BF<:BigFloat}(x::BF)
    isnan(x) && return x
    z = BF()
    ccall((:mpfr_sqrt, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
    if isnan(z)
        throw(DomainError())
    end
    return z
end

sqrt(x::BigInt) = sqrt(BF(x))

function ^{BF<:BigFloat}(x::BF, y::BF)
    z = BF()
    ccall((:mpfr_pow, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function ^{BF<:BigFloat}(x::BF, y::CulongMax)
    z = BF()
    ccall((:mpfr_pow_ui, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Culong, Int32), &z, &x, y, ROUNDING_MODE[])
    return z
end

function ^{BF<:BigFloat}(x::BF, y::ClongMax)
    z = BF()
    ccall((:mpfr_pow_si, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Clong, Int32), &z, &x, y, ROUNDING_MODE[])
    return z
end

function ^{BF<:BigFloat}(x::BF, y::BigInt)
    z = BF()
    ccall((:mpfr_pow_z, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BigInt}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

^{BF<:BigFloat}(x::BF, y::Integer)  = typemin(Clong)  <= y <= typemax(Clong)  ? x^Clong(y)  : x^BigInt(y)
^{BF<:BigFloat}(x::BF, y::Unsigned) = typemin(Culong) <= y <= typemax(Culong) ? x^Culong(y) : x^BigInt(y)

for f in (:exp, :exp2, :exp10, :expm1, :digamma, :erf, :erfc, :zeta,
          :cosh,:sinh,:tanh,:sech,:csch,:coth, :cbrt)
    @eval function $f{BF<:BigFloat}(x::BF)
        z = BF()
        ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
        return z
    end
end

# return log(2)
function big_ln2()
    c = BF()
    ccall((:mpfr_const_log2, :libmpfr), Cint, (Ptr{BF}, Int32),
          &c, MPFR.ROUNDING_MODE[])
    return c
end

function eta{BF<:BigFloat}(x::BF)
    x == 1 && return big_ln2()
    return -zeta(x) * expm1(big_ln2()*(1-x))
end

function airyai{BF<:BigFloat}(x::BF)
    z = BF()
    ccall((:mpfr_ai, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end
airy{BF<:BigFloat}(x::BF) = airyai(x)

function ldexp{BF<:BigFloat}(x::BF, n::Clong)
    z = BF()
    ccall((:mpfr_mul_2si, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Clong, Int32), &z, &x, n, ROUNDING_MODE[])
    return z
end
function ldexp{BF<:BigFloat}(x::BF, n::Culong)
    z = BF()
    ccall((:mpfr_mul_2ui, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Culong, Int32), &z, &x, n, ROUNDING_MODE[])
    return z
end
ldexp{BF<:BigFloat}(x::BF, n::ClongMax) = ldexp(x, convert(Clong, n))
ldexp{BF<:BigFloat}(x::BF, n::CulongMax) = ldexp(x, convert(Culong, n))
ldexp{BF<:BigFloat}(x::BF, n::Integer) = x*exp2(BF(n))

function besselj0{BF<:BigFloat}(x::BF)
    z = BF()
    ccall((:mpfr_j0, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function besselj1{BF<:BigFloat}(x::BF)
    z = BF()
    ccall((:mpfr_j1, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function besselj{BF<:BigFloat}(n::Integer, x::BF)
    z = BF()
    ccall((:mpfr_jn, :libmpfr), Int32, (Ptr{BF}, Clong, Ptr{BF}, Int32), &z, n, &x, ROUNDING_MODE[])
    return z
end

function bessely0{BF<:BigFloat}(x::BF)
    if x < 0
        throw(DomainError())
    end
    z = BF()
    ccall((:mpfr_y0, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function bessely1{BF<:BigFloat}(x::BF)
    if x < 0
        throw(DomainError())
    end
    z = BF()
    ccall((:mpfr_y1, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function bessely{BF<:BigFloat}(n::Integer, x::BF)
    if x < 0
        throw(DomainError())
    end
    z = BF()
    ccall((:mpfr_yn, :libmpfr), Int32, (Ptr{BF}, Clong, Ptr{BF}, Int32), &z, n, &x, ROUNDING_MODE[])
    return z
end

function factorial{BF<:BigFloat}(x::BF)
    if x < 0 || !isinteger(x)
        throw(DomainError())
    end
    ui = convert(Culong, x)
    z = BF()
    ccall((:mpfr_fac_ui, :libmpfr), Int32, (Ptr{BF}, Culong, Int32), &z, ui, ROUNDING_MODE[])
    return z
end

function hypot{BF<:BigFloat}(x::BF, y::BF)
    z = BF()
    ccall((:mpfr_hypot, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

for f in (:log, :log2, :log10)
    @eval function $f{BF<:BigFloat}(x::BF)
        if x < 0
            throw(DomainError())
        end
        z = BF()
        ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
        return z
    end
end

function log1p{BF<:BigFloat}(x::BF)
    if x < -1
        throw(DomainError())
    end
    z = BF()
    ccall((:mpfr_log1p, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
    return z
end

function max{BF<:BigFloat}(x::BF, y::BF)
    z = BF()
    ccall((:mpfr_max, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function min{BF<:BigFloat}(x::BF, y::BF)
    z = BF()
    ccall((:mpfr_min, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function modf{BF<:BigFloat}(x::BF)
    if isinf(x)
        return (BF(NaN), x)
    end
    zint = BF()
    zfloat = BF()
    ccall((:mpfr_modf, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &zint, &zfloat, &x, ROUNDING_MODE[])
    return (zfloat, zint)
end

function rem{BF<:BigFloat}(x::BF, y::BF)
    z = BF()
    ccall((:mpfr_fmod, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function sum{BF<:BigFloat}(arr::AbstractArray{BF})
    z = BF(0)
    for i in arr
        ccall((:mpfr_add, :libmpfr), Int32,
            (Ptr{BF}, Ptr{BF}, Ptr{BF}, Cint),
            &z, &z, &i, 0)
    end
    return z
end

# Functions for which NaN results are converted to DomainError, following Base
for f in (:sin,:cos,:tan,:sec,:csc,
          :acos,:asin,:atan,:acosh,:asinh,:atanh, :gamma)
    @eval begin
        function ($f){BF<:BigFloat}(x::BF)
            if isnan(x)
                return x
            end
            z = BF()
            ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32), &z, &x, ROUNDING_MODE[])
            if isnan(z)
                throw(DomainError())
            end
            return z
        end
    end
end

# log of absolute value of gamma function
const lgamma_signp = Array{Cint}(1)
function lgamma{BF<:BigFloat}(x::BF)
    z = BF()
    ccall((:mpfr_lgamma,:libmpfr), Cint, (Ptr{BF}, Ptr{Cint}, Ptr{BF}, Int32), &z, lgamma_signp, &x, ROUNDING_MODE[])
    return z
end

lgamma_r{BF<:BigFloat}(x::BF) = (lgamma(x), lgamma_signp[1])

function atan2{BF<:BigFloat}(y::BF, x::BF)
    z = BF()
    ccall((:mpfr_atan2, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &y, &x, ROUNDING_MODE[])
    return z
end

# Utility functions
=={BF<:BigFloat}(x::BF, y::BF) = ccall((:mpfr_equal_p, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &x, &y) != 0
<={BF<:BigFloat}(x::BF, y::BF) = ccall((:mpfr_lessequal_p, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &x, &y) != 0
>={BF<:BigFloat}(x::BF, y::BF) = ccall((:mpfr_greaterequal_p, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &x, &y) != 0
<{BF<:BigFloat}(x::BF, y::BF) = ccall((:mpfr_less_p, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &x, &y) != 0
>{BF<:BigFloat}(x::BF, y::BF) = ccall((:mpfr_greater_p, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &x, &y) != 0

function cmp{BF<:BigFloat}(x::BF, y::BigInt)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_z, :libmpfr), Int32, (Ptr{BF}, Ptr{BigInt}), &x, &y)
end
function cmp{BF<:BigFloat}(x::BF, y::ClongMax)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_si, :libmpfr), Int32, (Ptr{BF}, Clong), &x, y)
end
function cmp{BF<:BigFloat}(x::BF, y::CulongMax)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_ui, :libmpfr), Int32, (Ptr{BF}, Culong), &x, y)
end
cmp{BF<:BigFloat}(x::BF, y::Integer) = cmp(x,big(y))
cmp{BF<:BigFloat}(x::Integer, y::BF) = -cmp(y,x)

function cmp{BF<:BigFloat}(x::BF, y::CdoubleMax)
    (isnan(x) || isnan(y)) && throw(DomainError())
    ccall((:mpfr_cmp_d, :libmpfr), Int32, (Ptr{BF}, Cdouble), &x, y)
end
cmp{BF<:BigFloat}(x::CdoubleMax, y::BF) = -cmp(y,x)

=={BF<:BigFloat}(x::BF, y::Integer)   = !isnan(x) && cmp(x,y) == 0
=={BF<:BigFloat}(x::Integer, y::BF)   = y == x
=={BF<:BigFloat}(x::BF, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) == 0
=={BF<:BigFloat}(x::CdoubleMax, y::BF) = y == x

<{BF<:BigFloat}(x::BF, y::Integer)   = !isnan(x) && cmp(x,y) < 0
<{BF<:BigFloat}(x::Integer, y::BF)   = !isnan(y) && cmp(y,x) > 0
<{BF<:BigFloat}(x::BF, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) < 0
<{BF<:BigFloat}(x::CdoubleMax, y::BF) = !isnan(x) && !isnan(y) && cmp(y,x) > 0

<={BF<:BigFloat}(x::BF, y::Integer)   = !isnan(x) && cmp(x,y) <= 0
<={BF<:BigFloat}(x::Integer, y::BF)   = !isnan(y) && cmp(y,x) >= 0
<={BF<:BigFloat}(x::BF, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) <= 0
<={BF<:BigFloat}(x::CdoubleMax, y::BF) = !isnan(x) && !isnan(y) && cmp(y,x) >= 0

signbit{BF<:BigFloat}(x::BF) = ccall((:mpfr_signbit, :libmpfr), Int32, (Ptr{BF},), &x) != 0

function precision{BF<:BigFloat}(x::BF)  # precision of an object of type BF
    return ccall((:mpfr_get_prec, :libmpfr), Clong, (Ptr{BF},), &x)
end

precision{BF<:BigFloat}(::Type{BF}) = DEFAULT_PRECISION[end]  # precision of the type BF itself

"""
    setprecision([T=BF,] precision::Int)

Set the precision (in bits) to be used for `T` arithmetic.
"""
function setprecision(::Type{BigFloat}, precision::Int)
    if precision < 2
        throw(DomainError())
    end
    DEFAULT_PRECISION[end] = precision
end

setprecision(precision::Int) = setprecision(BigFloat, precision)

maxintfloat{BF<:BigFloat}(x::BF) = BF(2)^precision(x)
maxintfloat{BF<:BigFloat}(::Type{BF}) = BF(2)^precision(BF)

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

rounding_raw{BF<:BigFloat}(::Type{BF}) = ROUNDING_MODE[]
setrounding_raw{BF<:BigFloat}(::Type{BF},i::Integer) = ROUNDING_MODE[] = i

rounding{BF<:BigFloat}(::Type{BF}) = from_mpfr(rounding_raw(BF))
setrounding{BF<:BigFloat}(::Type{BF},r::RoundingMode) = setrounding_raw(BF,to_mpfr(r))

function copysign{BF<:BigFloat}(x::BF, y::BF)
    z = BF()
    ccall((:mpfr_copysign, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Ptr{BF}, Int32), &z, &x, &y, ROUNDING_MODE[])
    return z
end

function exponent{BF<:BigFloat}(x::BF)
    if x == 0 || !isfinite(x)
        throw(DomainError())
    end
    # The '- 1' is to make it work as Base.exponent
    return ccall((:mpfr_get_exp, :libmpfr), Clong, (Ptr{BF},), &x) - 1
end

function frexp{BF<:BigFloat}(x::BF)
    z = BF()
    c = Ref{Clong}()
    ccall((:mpfr_frexp, :libmpfr), Int32, (Ptr{Clong}, Ptr{BF}, Ptr{BF}, Cint), c, &z, &x, ROUNDING_MODE[])
    return (z, c[])
end

function significand{BF<:BigFloat}(x::BF)
    z = BF()
    c = Ref{Clong}()
    ccall((:mpfr_frexp, :libmpfr), Int32, (Ptr{Clong}, Ptr{BF}, Ptr{BF}, Cint), c, &z, &x, ROUNDING_MODE[])
    # Double the significand to make it work as Base.significand
    ccall((:mpfr_mul_si, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Clong, Int32), &z, &z, 2, ROUNDING_MODE[])
    return z
end

function isinteger{BF<:BigFloat}(x::BF)
    return ccall((:mpfr_integer_p, :libmpfr), Int32, (Ptr{BF},), &x) != 0
end

for f in (:ceil, :floor, :trunc)
    @eval begin
        function ($f){BF<:BigFloat}(x::BF)
            z = BF()
            ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &x)
            return z
        end
    end
end

function round{BF<:BigFloat}(x::BF)
    z = BF()
    ccall((:mpfr_rint, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Cint), &z, &x, ROUNDING_MODE[])
    return z
end
function round{BF<:BigFloat}(x::BF,::RoundingMode{:NearestTiesAway})
    z = BF()
    ccall((:mpfr_round, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}), &z, &x)
    return z
end

function isinf{BF<:BigFloat}(x::BF)
    return ccall((:mpfr_inf_p, :libmpfr), Int32, (Ptr{BF},), &x) != 0
end

function isnan{BF<:BigFloat}(x::BF)
    return ccall((:mpfr_nan_p, :libmpfr), Int32, (Ptr{BF},), &x) != 0
end

isfinite{BF<:BigFloat}(x::BF) = !isinf(x) && !isnan(x)

@eval typemax(::Type{BigFloat}) = $(BigFloat( Inf))
@eval typemin(::Type{BigFloat}) = $(BigFloat(-Inf))
typemax{BF<:BigFloat}(::Type{BF}) = BF( Inf)
typemin{BF<:BigFloat}(::Type{BF}) = BF(-Inf)

function nextfloat{BF<:BigFloat}(x::BF)
    z = BF()
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32),
          &z, &x, ROUNDING_MODE[])
    ccall((:mpfr_nextabove, :libmpfr), Int32, (Ptr{BF},), &z) != 0
    return z
end

function prevfloat{BF<:BigFloat}(x::BF)
    z = BF()
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BF}, Ptr{BF}, Int32),
          &z, &x, ROUNDING_MODE[])
    ccall((:mpfr_nextbelow, :libmpfr), Int32, (Ptr{BF},), &z) != 0
    return z
end

eps{BF<:BigFloat}(::Type{BF}) = nextfloat(BF(1)) - BF(1)

realmin{BF<:BigFloat}(::Type{BF}) = nextfloat(zero(BF))
realmax{BF<:BigFloat}(::Type{BF}) = prevfloat(BF(Inf))

"""
    setprecision(f::Function, [T=BF,] precision::Integer)

Change the `T` arithmetic precision (in bits) for the duration of `f`.
It is logically equivalent to:

    old = precision(BF)
    setprecision(BF, precision)
    f()
    setprecision(BF, old)

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

function string(x::BigFloat)
    # In general, the number of decimal places needed to read back the number exactly
    # is, excluding the most significant, ceil(log(10, 2^precision(x)))
    k = ceil(Int32, precision(x) * 0.3010299956639812)
    lng = k + Int32(8) # Add space for the sign, the most significand digit, the dot and the exponent
    buf = Array{UInt8}(lng + 1)
    # format strings are guaranteed to contain no NUL, so we don't use Cstring
    lng = ccall((:mpfr_snprintf,:libmpfr), Int32, (Ptr{UInt8}, Culong, Ptr{UInt8}, Ptr{BigFloat}...), buf, lng + 1, "%.Re", &x)
    if lng < k + 5 # print at least k decimal places
        lng = ccall((:mpfr_sprintf,:libmpfr), Int32, (Ptr{UInt8}, Ptr{UInt8}, Ptr{BigFloat}...), buf, "%.$(k)Re", &x)
    elseif lng > k + 8
        buf = Array{UInt8}(lng + 1)
        lng = ccall((:mpfr_snprintf,:libmpfr), Int32, (Ptr{UInt8}, Culong, Ptr{UInt8}, Ptr{BigFloat}...), buf, lng + 1, "%.Re", &x)
    end
    n = (1 <= x < 10 || -10 < x <= -1 || x == 0) ? lng - 4 : lng
    return String(buf[1:n])
end

function string{BF<:BigFloat}(x::BF)
    # In general, the number of decimal places needed to read back the number exactly
    # is, excluding the most significant, ceil(log(10, 2^precision(x)))
    k = ceil(Int32, precision(x) * 0.3010299956639812)
    lng = k + Int32(8) # Add space for the sign, the most significand digit, the dot and the exponent
    buf = Array{UInt8}(lng + 1)
    # format strings are guaranteed to contain no NUL, so we don't use Cstring
    lng = ccall((:mpfr_snprintf,:libmpfr), Int32, (Ptr{UInt8}, Culong, Ptr{UInt8}, Ptr{BF}...), buf, lng + 1, "%.Re", &x)
    if lng < k + 5 # print at least k decimal places
        lng = ccall((:mpfr_sprintf,:libmpfr), Int32, (Ptr{UInt8}, Ptr{UInt8}, Ptr{BF}...), buf, "%.$(k)Re", &x)
    elseif lng > k + 8
        buf = Array{UInt8}(lng + 1)
        lng = ccall((:mpfr_snprintf,:libmpfr), Int32, (Ptr{UInt8}, Culong, Ptr{UInt8}, Ptr{BF}...), buf, lng + 1, "%.Re", &x)
    end
    n = (1 <= x < 10 || -10 < x <= -1 || x == 0) ? lng - 4 : lng
    return String(buf[1:n])
end

print(io::IO, b::BigFloat) = print(io, string(b))
show(io::IO, b::BigFloat) = print(io, string(b))
print{BF<:BigFloat}(io::IO, b::BF) = print(io, string(b))
show{BF<:BigFloat}(io::IO, b::BF) = print(io, string(b))


# get/set exponent min/max
get_emax() = ccall((:mpfr_get_emax, :libmpfr), Clong, ())
get_emax_min() = ccall((:mpfr_get_emax_min, :libmpfr), Clong, ())
get_emax_max() = ccall((:mpfr_get_emax_max, :libmpfr), Clong, ())

get_emin() = ccall((:mpfr_get_emin, :libmpfr), Clong, ())
get_emin_min() = ccall((:mpfr_get_emin_min, :libmpfr), Clong, ())
get_emin_max() = ccall((:mpfr_get_emin_max, :libmpfr), Clong, ())

set_emax!(x) = ccall((:mpfr_set_emax, :libmpfr), Void, (Clong,), x)
set_emin!(x) = ccall((:mpfr_set_emin, :libmpfr), Void, (Clong,), x)

function Base.deepcopy_internal(x::BigFloat, stackdict::ObjectIdDict)
    if haskey(stackdict, x)
        return stackdict[x]
    end
    N = precision(x)
    y = BigFloat(zero(Clong), zero(Cint), zero(Clong), C_NULL)
    ccall((:mpfr_init2,:libmpfr), Void, (Ptr{BigFloat}, Clong), &y, N)
    finalizer(y, cglobal((:mpfr_clear, :libmpfr)))
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Int32),
          &y, &x, ROUNDING_MODE[])
    stackdict[x] = y
    return y
end

convert{BF<:BigFloat}(::Type{BF}, x::BF) = x
function convert{BF1<:BigFloat, BF2<:BigFloat}(::Type{BF1}, x::BF2)
    return parse(BF1, string(x))
end
function promote_rule{BF1<:BigFloat, BF2<:BigFloat}(::Type{BF1}, ::Type{BF2})
    return ifelse(precision(BF1) < precision(BF2), BF2, BF1)
end

(+){BF1<:BigFloat, BF2<:BigFloat}(a::BF1, b::BF2) = (+)(promote(a,b)...)
(-){BF1<:BigFloat, BF2<:BigFloat}(a::BF1, b::BF2) = (-)(promote(a,b)...)
(*){BF1<:BigFloat, BF2<:BigFloat}(a::BF1, b::BF2) = (*)(promote(a,b)...)
(/){BF1<:BigFloat, BF2<:BigFloat}(a::BF1, b::BF2) = (/)(promote(a,b)...)
(hypot){BF1<:BigFloat, BF2<:BigFloat}(a::BF1, b::BF2) = (hypot)(promote(a,b)...)
(^){BF1<:BigFloat, BF2<:BigFloat}(a::BF1, b::BF2) = (^)(promote(a,b)...)


end #module

