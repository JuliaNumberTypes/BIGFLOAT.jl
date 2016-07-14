
# total ordering
=={P}(x::BigFP{P}, y::BigFP{P}) = ccall((:mpfr_equal_p, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &x, &y) != 0
!={P}(x::BigFP{P}, y::BigFP{P}) = ccall((:mpfr_equal_p, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &x, &y) == 0
<={P}(x::BigFP{P}, y::BigFP{P}) = ccall((:mpfr_lessequal_p, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &x, &y) != 0
>={P}(x::BigFP{P}, y::BigFP{P}) = ccall((:mpfr_greaterequal_p, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &x, &y) != 0
 <{P}(x::BigFP{P}, y::BigFP{P}) = ccall((:mpfr_less_p, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &x, &y) != 0
 >{P}(x::BigFP{P}, y::BigFP{P}) = ccall((:mpfr_greater_p, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &x, &y) != 0

# utility functions for comparisons
function cmp{P}(x::BigFP{P}, y::BigInt)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_z, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigInt}), &x, &y)
end
function cmp{P}(x::BigFP{P}, y::ClongMax)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_si, :libmpfr), Int32, (Ptr{BigFP{P}}, Clong), &x, y)
end
function cmp{P}(x::BigFP{P}, y::CulongMax)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_ui, :libmpfr), Int32, (Ptr{BigFP{P}}, Culong), &x, y)
end
cmp{P}(x::BigFP{P}, y::Integer) = cmp(x,big(y))
cmp{P}(x::Integer, y::BigFP{P}) = -cmp(y,x)

function cmp{P}(x::BigFP{P}, y::CdoubleMax)
    (isnan(x) || isnan(y)) && throw(DomainError())
    ccall((:mpfr_cmp_d, :libmpfr), Int32, (Ptr{BigFP{P}}, Cdouble), &x, y)
end
cmp{P}(x::CdoubleMax, y::BigFP{P}) = -cmp(y,x)

# intertype compare

=={P}(x::BigFP{P}, y::Integer)    = !isnan(x) && cmp(x,y) == 0
=={P}(x::Integer, y::BigFP{P})    = y == x
=={P}(x::BigFP{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) == 0
=={P}(x::CdoubleMax, y::BigFP{P}) = y == x

!={P}(x::BigFP{P}, y::Integer)    = isnan(x) || cmp(x,y) != 0
!={P}(x::Integer, y::BigFP{P})    = y != x
!={P}(x::BigFP{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) != 0
!={P}(x::CdoubleMax, y::BigFP{P}) = y != x

 <{P}(x::BigFP{P}, y::Integer)    = !isnan(x) && cmp(x,y) < 0
 <{P}(x::Integer, y::BigFP{P})    = !isnan(y) && cmp(y,x) > 0
 <{P}(x::BigFP{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) < 0
 <{P}(x::CdoubleMax, y::BigFP{P}) = !isnan(x) && !isnan(y) && cmp(y,x) > 0

<={P}(x::BigFP{P}, y::Integer)    = !isnan(x) && cmp(x,y) <= 0
<={P}(x::Integer, y::BigFP{P})    = !isnan(y) && cmp(y,x) >= 0
<={P}(x::BigFP{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) <= 0
<={P}(x::CdoubleMax, y::BigFP{P}) = !isnan(x) && !isnan(y) && cmp(y,x) >= 0

 >{P}(x::BigFP{P}, y::Integer)    = !isnan(x) && cmp(x,y) > 0
 >{P}(x::Integer, y::BigFP{P})    = !isnan(y) && cmp(y,x) < 0
 >{P}(x::BigFP{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) > 0
 >{P}(x::CdoubleMax, y::BigFP{P}) = !isnan(x) && !isnan(y) && cmp(y,x) < 0

>={P}(x::BigFP{P}, y::Integer)    = !isnan(x) && cmp(x,y) >= 0
>={P}(x::Integer, y::BigFP{P})    = !isnan(y) && cmp(y,x) <= 0
>={P}(x::BigFP{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) >= 0
>={P}(x::CdoubleMax, y::BigFP{P}) = !isnan(x) && !isnan(y) && cmp(y,x) <= 0
