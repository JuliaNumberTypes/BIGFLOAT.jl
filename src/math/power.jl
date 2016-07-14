function (^){P}(x::BigFloat{P}, y::BigFloat{P})
    z = BigFloat{P}()
    ccall((:mpfr_pow, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end

function (^){P}(x::BigFloat{P}, y::CulongMax)
    z = BigFloat{P}()
    ccall((:mpfr_pow_ui, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Culong, Int32), &z, &x, y, ROUNDING_MODE[end])
    return z
end

function (^){P}(x::BigFloat{P}, y::ClongMax)
    z = BigFloat{P}()
    ccall((:mpfr_pow_si, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Clong, Int32), &z, &x, y, ROUNDING_MODE[end])
    return z
end

function (^){P}(x::BigFloat{P}, y::BigInt)
    z = BigFloat{P}()
    ccall((:mpfr_pow_z, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Ptr{BigInt}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end

(^){P}(x::BigFloat{P}, y::Integer)  = typemin(Clong)  <= y <= typemax(Clong)  ? x^Clong(y)  : x^BigInt(y)
(^){P}(x::BigFloat{P}, y::Unsigned) = typemin(Culong) <= y <= typemax(Culong) ? x^Culong(y) : x^BigInt(y)
