function (^){P}(x::BigFP{P}, y::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_pow, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end

function (^){P}(x::BigFP{P}, y::CulongMax)
    z = BigFP{P}()
    ccall((:mpfr_pow_ui, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Culong, Int32), &z, &x, y, ROUNDING_MODE[end])
    return z
end

function (^){P}(x::BigFP{P}, y::ClongMax)
    z = BigFP{P}()
    ccall((:mpfr_pow_si, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Clong, Int32), &z, &x, y, ROUNDING_MODE[end])
    return z
end

function (^){P}(x::BigFP{P}, y::BigInt)
    z = BigFP{P}()
    ccall((:mpfr_pow_z, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigInt}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end

(^){P}(x::BigFP{P}, y::Integer)  = typemin(Clong)  <= y <= typemax(Clong)  ? x^Clong(y)  : x^BigInt(y)
(^){P}(x::BigFP{P}, y::Unsigned) = typemin(Culong) <= y <= typemax(Culong) ? x^Culong(y) : x^BigInt(y)
