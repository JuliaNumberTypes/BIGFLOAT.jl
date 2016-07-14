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

rounding_raw(::Type{BigFP}) = ROUNDING_MODE[end]
setrounding_raw(::Type{BigFP},i::Integer) = ROUNDING_MODE[end] = i
rounding_raw{P}(::Type{BigFP{P}}) = rounding_raw(BigFP)
setrounding_raw{P}(::Type{BigFP{P}},i::Integer) = setrounding_raw(BigFP)

rounding(::Type{BigFP}) = from_mpfr(rounding_raw(BigFP))
setrounding(::Type{BigFP},r::RoundingMode) = setrounding_raw(BigFP,to_mpfr(r))
rounding{P}(::Type{BigFP{P}}) = rounding(BigFP)
setrounding{P}(::Type{BigFP{P}},r::RoundingMode) = setrounding(BigFP,to_mpfr(r))


function round{P}(x::BigFP{P})
    z = BigFP{P}()
    ccall((:mpfr_rint, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}, Cint), &z, &x, ROUNDING_MODE[end])
    return z
end
function round{P}(x::BigFP{P},::RoundingMode{:NearestTiesAway})
    z = BigFP{P}()
    ccall((:mpfr_round, :libmpfr), Int32, (Ptr{BigFP{P}}, Ptr{BigFP{P}}), &z, &x)
    return z
end
