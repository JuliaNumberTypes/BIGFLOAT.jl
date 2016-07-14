
for f in (:exp, :exp2, :exp10, :expm1, :digamma, :erf, :erfc, :zeta,
          :cosh,:sinh,:tanh,:sech,:csch,:coth, :cbrt)
    @eval function ($f){P}(x::BigFP{P})
        z = BigFP{P}()
        ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
        return z
    end
end

# return log(2)
function big_ln2()
    c = BigFP{P}()
    ccall((:mpfr_const_log2, :libmpfr), Cint, (Ptr{BigFP{P}}, Int32),
          &c, MPFR.ROUNDING_MODE[end])
    return c
end

function eta{P}(x::BigFP{P})
    x == 1 && return big_ln2()
    return -zeta(x) * expm1(big_ln2()*(1-x))
end

function airyai{P}(x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_ai, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
    return z
end
airy(x::BigFP{P}) = airyai(x)

function ldexp{P}(x::BigFP{P}, n::Clong)
    z = BigFP{P}()
    ccall((:mpfr_mul_2si, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Clong, Int32), &z, &x, n, ROUNDING_MODE[end])
    return z
end
function ldexp{P}(x::BigFP{P}, n::Culong)
    z = BigFP{P}()
    ccall((:mpfr_mul_2ui, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Culong, Int32), &z, &x, n, ROUNDING_MODE[end])
    return z
end
ldexp{P}(x::BigFP{P}, n::ClongMax) = ldexp(x, convert(Clong, n))
ldexp{P}(x::BigFP{P}, n::CulongMax) = ldexp(x, convert(Culong, n))
ldexp{P}(x::BigFP{P}, n::Integer) = x*exp2(BigFP{P}(n))

function besselj0{P}(x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_j0, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
    return z
end

function besselj1{P}(x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_j1, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
    return z
end

function besselj{P}(n::Integer, x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_jn, :libmpfr), Int32, (Ptr{BigFP{P}}, Clong, Ptr{BigFP{P}}, Int32), &z, n, &x, ROUNDING_MODE[end])
    return z
end

function bessely0{P}(x::BigFP{P})
    if x < 0
        throw(DomainError())
    end
    z = BigFP{P}()
    ccall((:mpfr_y0, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
    return z
end

function bessely1{P}(x::BigFP{P})
    if x < 0
        throw(DomainError())
    end
    z = BigFP{P}()
    ccall((:mpfr_y1, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
    return z
end

function bessely{P}(n::Integer, x::BigFP{P})
    if x < 0
        throw(DomainError())
    end
    z = BigFP{P}()
    ccall((:mpfr_yn, :libmpfr), Int32, (Ptr{BigFP{P}}, Clong, Ptr{BigFP{P}}, Int32), &z, n, &x, ROUNDING_MODE[end])
    return z
end

function factorial{P}(x::BigFP{P})
    if x < 0 || !isinteger(x)
        throw(DomainError())
    end
    ui = convert(Culong, x)
    z = BigFP{P}()
    ccall((:mpfr_fac_ui, :libmpfr), Int32, (Ptr{BigFP{P}}, Culong, Int32), &z, ui, ROUNDING_MODE[end])
    return z
end

for f in (:log, :log2, :log10)
    @eval function ($f){P}(x::BigFP{P})
        if x < 0
            throw(DomainError())
        end
        z = BigFP{P}()
        ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
        return z
    end
end

function log1p{P}(x::BigFP{P})
    if x < -1
        throw(DomainError())
    end
    z = BigFP{P}()
    ccall((:mpfr_log1p, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
    return z
end


function modf{P}(x::BigFP{P})
    if isinf(x)
        return (BigFP{P}(NaN), x)
    end
    zint = BigFP{P}()
    zfloat = BigFP{P}()
    ccall((:mpfr_modf, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &zint, &zfloat, &x, ROUNDING_MODE[end])
    return (zfloat, zint)
end

function rem{P}(x::BigFP{P}, y::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_fmod, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end


# Functions for which NaN results are converted to DomainError, following Base
for f in (:sin,:cos,:tan,:sec,:csc,
          :acos,:asin,:atan,:acosh,:asinh,:atanh, :gamma)
    @eval begin
        function ($f){P}(x::BigFP{P})
            if isnan(x)
                return x
            end
            z = BigFP{P}()
            ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
            if isnan(z)
                throw(DomainError())
            end
            return z
        end
    end
end

# log of absolute value of gamma function
const lgamma_signp = Array{Cint}(1)
function lgamma{P}(x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_lgamma,:libmpfr), Cint, (Ptr{BigFP{P}}, Ptr{Cint}, Ptr{BigFP{P}}, Int32), &z, lgamma_signp, &x, ROUNDING_MODE[end])
    return z
end

lgamma_r{P}(x::BigFP{P}) = (lgamma(x), lgamma_signp[1])

function atan2{P}(y::BigFP{P}, x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_atan2, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &y, &x, ROUNDING_MODE[end])
    return z
end
