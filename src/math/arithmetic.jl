# Basic arithmetic without promotion
for (fJ, fC) in ((:+,:add), (:*,:mul))
    @eval begin
        # BigFP{P}
        function ($fJ){P}(x::BigFP{P}, y::BigFP{P})
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
            return z
        end

        # Unsigned Integer
        function ($fJ){P}(x::BigFP{P}, c::CulongMax)
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC,:_ui)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Culong, Int32), &z, &x, c, ROUNDING_MODE[end])
            return z
        end
        ($fJ){P}(c::CulongMax, x::BigFP{P}) = ($fJ){P}(x,c)

        # Signed Integer
        function ($fJ){P}(x::BigFP{P}, c::ClongMax)
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC,:_si)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Clong, Int32), &z, &x, c, ROUNDING_MODE[end])
            return z
        end
        ($fJ){P}(c::ClongMax, x::BigFP{P}) = ($fJ){P}(x,c)

        # Float32/Float64
        function ($fJ){P}(x::BigFP{P}, c::CdoubleMax)
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC,:_d)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Cdouble, Int32), &z, &x, c, ROUNDING_MODE[end])
            return z
        end
        ($fJ){P}(c::CdoubleMax, x::BigFP{P}) = ($fJ){P}(x,c)

        # BigInt
        function ($fJ){P}(x::BigFP{P}, c::BigInt)
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC,:_z)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigInt}, Int32), &z, &x, &c, ROUNDING_MODE[end])
            return z
        end
        ($fJ){P}(c::BigInt, x::BigFP{P}) = ($fJ){P}(x,c)
    end
end

for (fJ, fC) in ((:-,:sub), (:/,:div))
    @eval begin
        # BigFP{P}
        function ($fJ){P}(x::BigFP{P}, y::BigFP{P})
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
            return z
        end

        # Unsigned Int
        function ($fJ){P}(x::BigFP{P}, c::CulongMax)
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC,:_ui)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Culong, Int32), &z, &x, c, ROUNDING_MODE[end])
            return z
        end
        function ($fJ){P}(c::CulongMax, x::BigFP{P})
            z = BigFP{P}()
            ccall(($(string(:mpfr_,:ui_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Culong, Ptr{BigFP{P}}, Int32), &z, c, &x, ROUNDING_MODE[end])
            return z
        end

        # Signed Integer
        function ($fJ){P}(x::BigFP{P}, c::ClongMax)
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC,:_si)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Clong, Int32), &z, &x, c, ROUNDING_MODE[end])
            return z
        end
        function ($fJ){P}(c::ClongMax, x::BigFP{P})
            z = BigFP{P}()
            ccall(($(string(:mpfr_,:si_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Clong, Ptr{BigFP{P}}, Int32), &z, c, &x, ROUNDING_MODE[end])
            return z
        end

        # Float32/Float64
        function ($fJ){P}(x::BigFP{P}, c::CdoubleMax)
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC,:_d)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Cdouble, Int32), &z, &x, c, ROUNDING_MODE[end])
            return z
        end
        function ($fJ){P}(c::CdoubleMax, x::BigFP{P})
            z = BigFP{P}()
            ccall(($(string(:mpfr_,:d_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Cdouble, Ptr{BigFP{P}}, Int32), &z, c, &x, ROUNDING_MODE[end])
            return z
        end

        # BigInt
        function ($fJ){P}(x::BigFP{P}, c::BigInt)
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC,:_z)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigInt}, Int32), &z, &x, &c, ROUNDING_MODE[end])
            return z
        end
        # no :mpfr_z_div function
    end
end

function (-){P}(c::BigInt, x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_z_sub, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigInt}, Ptr{BigFP{P}}, Int32), &z, &c, &x, ROUNDING_MODE[end])
    return z
end

function fma{P}(x::BigFP{P}, y::BigFP{P}, z::BigFP{P})
    r = BigFP{P}()
    ccall(("mpfr_fma",:libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &r, &x, &y, &z, ROUNDING_MODE[end])
    return r
end

# div
# BigFP{P}
function div{P}(x::BigFP{P}, y::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_div,:libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, &y, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &z, &z)
    return z
end

# Unsigned Int
function div{P}(x::BigFP{P}, c::CulongMax)
    z = BigFP{P}()
    ccall((:mpfr_div_ui, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Culong, Int32), &z, &x, c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &z, &z)
    return z
end
function div{P}(c::CulongMax, x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_ui_div, :libmpfr), Int32, (Ptr{BigFP{P}}, Culong, Ptr{BigFP{P}}, Int32), &z, c, &x, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &z, &z)
    return z
end

# Signed Integer
function div{P}(x::BigFP{P}, c::ClongMax)
    z = BigFP{P}()
    ccall((:mpfr_div_si, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Clong, Int32), &z, &x, c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &z, &z)
    return z
end
function div{P}(c::ClongMax, x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_si_div, :libmpfr), Int32, (Ptr{BigFP{P}}, Clong, Ptr{BigFP{P}}, Int32), &z, c, &x, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &z, &z)
    return z
end

# Float32/Float64
function div{P}(x::BigFP{P}, c::CdoubleMax)
    z = BigFP{P}()
    ccall((:mpfr_div_d, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Cdouble, Int32), &z, &x, c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &z, &z)
    return z
end
function div{P}(c::CdoubleMax, x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_d_div, :libmpfr), Int32, (Ptr{BigFP{P}}, Cdouble, Ptr{BigFP{P}}, Int32), &z, c, &x, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &z, &z)
    return z
end

# BigInt
function div{P}(x::BigFP{P}, c::BigInt)
    z = BigFP{P}()
    ccall((:mpfr_div_z, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigInt}, Int32), &z, &x, &c, to_mpfr(RoundToZero))
    ccall((:mpfr_trunc, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &z, &z)
    return z
end


# More efficient commutative operations
for (fJ, fC, fI) in ((:+, :add, 0), (:*, :mul, 1))
    @eval begin
        function ($fJ){P}(a::BigFP{P}, b::BigFP{P}, c::BigFP{P})
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &a, &b, ROUNDING_MODE[end])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &z, &c, ROUNDING_MODE[end])
            return z
        end
        function ($fJ){P}(a::BigFP{P}, b::BigFP{P}, c::BigFP{P}, d::BigFP{P})
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &a, &b, ROUNDING_MODE[end])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &z, &c, ROUNDING_MODE[end])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &z, &d, ROUNDING_MODE[end])
            return z
        end
        function ($fJ){P}(a::BigFP{P}, b::BigFP{P}, c::BigFP{P}, d::BigFP{P}, e::BigFP{P})
            z = BigFP{P}()
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &a, &b, ROUNDING_MODE[end])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &z, &c, ROUNDING_MODE[end])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &z, &d, ROUNDING_MODE[end])
            ccall(($(string(:mpfr_,fC)), :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &z, &e, ROUNDING_MODE[end])
            return z
        end
    end
end

function (-){P}(x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_neg, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
    return z
end

function sqrt{P}(x::BigFP{P})
    isnan(x) && return x
    z = BigFP{P}()
    ccall((:mpfr_sqrt, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, ROUNDING_MODE[end])
    if isnan(z)
        throw(DomainError())
    end
    return z
end

sqrt(x::BigInt) = sqrt(BigFP{P}(x))

function hypot{P}(x::BigFP{P}, y::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_hypot, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Int32), &z, &x, &y, ROUNDING_MODE[end])
    return z
end



function sum{P}(arr::AbstractArray{BigFP{P}})
    z = BigFP{P}(0)
    for i in arr
        ccall((:mpfr_add, :libmpfr), Int32,
            (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Ptr{BigFP{P}}, Cint),
            &z, &z, &i, 0)
    end
    return z
end
