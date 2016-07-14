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

rounding_raw(::Type{BigFloat}) = ROUNDING_MODE[end]
setrounding_raw(::Type{BigFloat},i::Integer) = ROUNDING_MODE[end] = i
rounding_raw{P}(::Type{BigFloat{P}}) = rounding_raw(BigFloat)
setrounding_raw{P}(::Type{BigFloat{P}},i::Integer) = setrounding_raw(BigFloat)

rounding(::Type{BigFloat}) = from_mpfr(rounding_raw(BigFloat))
setrounding(::Type{BigFloat},r::RoundingMode) = setrounding_raw(BigFloat,to_mpfr(r))
rounding{P}(::Type{BigFloat{P}}) = rounding(BigFloat)
setrounding{P}(::Type{BigFloat{P}},r::RoundingMode) = setrounding(BigFloat,to_mpfr(r))


function round{P}(x::BigFloat{P})
    z = BigFloat{P}()
    ccall((:mpfr_rint, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}, Cint), &z, &x, ROUNDING_MODE[end])
    return z
end
function round{P}(x::BigFloat{P},::RoundingMode{:NearestTiesAway})
    z = BigFloat{P}()
    ccall((:mpfr_round, :libmpfr), Int32, (Ptr{BigFloat{P}}, Ptr{BigFloat{P}}), &z, &x)
    return z
end
