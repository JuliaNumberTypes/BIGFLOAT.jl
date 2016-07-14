
# total ordering
=={P}(x::BigFloat{P}, y::BigFloat{P}) = ccall((:mpfr_equal_p, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}), &x, &y) != 0
!={P}(x::BigFloat{P}, y::BigFloat{P}) = ccall((:mpfr_equal_p, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}), &x, &y) == 0
<={P}(x::BigFloat{P}, y::BigFloat{P}) = ccall((:mpfr_lessequal_p, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}), &x, &y) != 0
>={P}(x::BigFloat{P}, y::BigFloat{P}) = ccall((:mpfr_greaterequal_p, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}), &x, &y) != 0
 <{P}(x::BigFloat{P}, y::BigFloat{P}) = ccall((:mpfr_less_p, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}), &x, &y) != 0
 >{P}(x::BigFloat{P}, y::BigFloat{P}) = ccall((:mpfr_greater_p, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}), &x, &y) != 0

# utility functions for comparisons
function cmp{P}(x::BigFloat{P}, y::BigInt)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_z, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigInt}), &x, &y)
end
function cmp{P}(x::BigFloat{P}, y::ClongMax)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_si, :libmpfr), Int32, (Ptr{BigFloat{P}}, Clong), &x, y)
end
function cmp{P}(x::BigFloat{P}, y::CulongMax)
    isnan(x) && throw(DomainError())
    ccall((:mpfr_cmp_ui, :libmpfr), Int32, (Ptr{BigFloat{P}}, Culong), &x, y)
end
cmp{P}(x::BigFloat{P}, y::Integer) = cmp(x,big(y))
cmp{P}(x::Integer, y::BigFloat{P}) = -cmp(y,x)

function cmp{P}(x::BigFloat{P}, y::CdoubleMax)
    (isnan(x) || isnan(y)) && throw(DomainError())
    ccall((:mpfr_cmp_d, :libmpfr), Int32, (Ptr{BigFloat{P}}, Cdouble), &x, y)
end
cmp{P}(x::CdoubleMax, y::BigFloat{P}) = -cmp(y,x)

# intertype compare

=={P}(x::BigFloat{P}, y::Integer)    = !isnan(x) && cmp(x,y) == 0
=={P}(x::Integer, y::BigFloat{P})    = y == x
=={P}(x::BigFloat{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) == 0
=={P}(x::CdoubleMax, y::BigFloat{P}) = y == x

!={P}(x::BigFloat{P}, y::Integer)    = isnan(x) || cmp(x,y) != 0
!={P}(x::Integer, y::BigFloat{P})    = y != x
!={P}(x::BigFloat{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) != 0
!={P}(x::CdoubleMax, y::BigFloat{P}) = y != x

 <{P}(x::BigFloat{P}, y::Integer)    = !isnan(x) && cmp(x,y) < 0
 <{P}(x::Integer, y::BigFloat{P})    = !isnan(y) && cmp(y,x) > 0
 <{P}(x::BigFloat{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) < 0
 <{P}(x::CdoubleMax, y::BigFloat{P}) = !isnan(x) && !isnan(y) && cmp(y,x) > 0

<={P}(x::BigFloat{P}, y::Integer)    = !isnan(x) && cmp(x,y) <= 0
<={P}(x::Integer, y::BigFloat{P})    = !isnan(y) && cmp(y,x) >= 0
<={P}(x::BigFloat{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) <= 0
<={P}(x::CdoubleMax, y::BigFloat{P}) = !isnan(x) && !isnan(y) && cmp(y,x) >= 0

 >{P}(x::BigFloat{P}, y::Integer)    = !isnan(x) && cmp(x,y) > 0
 >{P}(x::Integer, y::BigFloat{P})    = !isnan(y) && cmp(y,x) < 0
 >{P}(x::BigFloat{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) > 0
 >{P}(x::CdoubleMax, y::BigFloat{P}) = !isnan(x) && !isnan(y) && cmp(y,x) < 0

>={P}(x::BigFloat{P}, y::Integer)    = !isnan(x) && cmp(x,y) >= 0
>={P}(x::Integer, y::BigFloat{P})    = !isnan(y) && cmp(y,x) <= 0
>={P}(x::BigFloat{P}, y::CdoubleMax) = !isnan(x) && !isnan(y) && cmp(x,y) >= 0
>={P}(x::CdoubleMax, y::BigFloat{P}) = !isnan(x) && !isnan(y) && cmp(y,x) <= 0
