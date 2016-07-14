
convert(::Type{BigFP}, x::BigFP) = x

# convert to BigFP
for (fJ, fC) in ((:si,:Clong), (:ui,:Culong), (:d,:Float64))
    @eval begin
        function convert(::Type{BigFP}, x::($fC))
            z = BigFP()
            ccall(($(string(:mpfr_set_,fJ)), :libmpfr), Int32, (Ptr{BigFP}, ($fC), Int32), &z, x, ROUNDING_MODE[end])
            return z
        end
    end
end

function convert(::Type{BigFP}, x::BigInt)
    z = BigFP()
    ccall((:mpfr_set_z, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigInt}, Int32), &z, &x, ROUNDING_MODE[end])
    return z
end

convert(::Type{BigFP}, x::Integer) = BigFP(BigInt(x))

convert(::Type{BigFP}, x::Union{Bool,Int8,Int16,Int32}) = BigFP(convert(Clong,x))
convert(::Type{BigFP}, x::Union{UInt8,UInt16,UInt32}) = BigFP(convert(Culong,x))

convert(::Type{BigFP}, x::Union{Float16,Float32}) = BigFP(Float64(x))
convert(::Type{BigFP}, x::Rational) = BigFP(num(x)) / BigFP(den(x))

function tryparse(::Type{BigFP}, s::AbstractString, base::Int=0)
    z = BigFP()
    err = ccall((:mpfr_set_str, :libmpfr), Int32, (Ptr{BigFP}, Cstring, Int32, Int32), &z, s, base, ROUNDING_MODE[end])
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
    ccall((:mpfr_get_d,:libmpfr), Float64, (Ptr{BigFP},Int32), &x, ROUNDING_MODE[end])
convert(::Type{Float32}, x::BigFP) =
    ccall((:mpfr_get_flt,:libmpfr), Float32, (Ptr{BigFP},Int32), &x, ROUNDING_MODE[end])
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
