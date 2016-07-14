
convert(::Type{BigFloat}, x::BigFloat) = x

# convert to BigFloat
for (fJ, fC) in ((:si,:Clong), (:ui,:Culong), (:d,:Float64))
    @eval begin
        function convert(::Type{BigFloat}, x::($fC))
            z = BigFloat()
            ccall(($(string(:mpfr_set_,fJ)), :libmpfr), Int32, (Ptr{BigFloat}, ($fC), Int32), &z, x, ROUNDING_MODE[end])
            return z
        end
    end
end

function convert(::Type{BigFloat}, x::BigInt)
    z = BigFloat()
    ccall((:mpfr_set_z, :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigInt}, Int32), &z, &x, ROUNDING_MODE[end])
    return z
end

convert(::Type{BigFloat}, x::Integer) = BigFloat(BigInt(x))

convert(::Type{BigFloat}, x::Union{Bool,Int8,Int16,Int32}) = BigFloat(convert(Clong,x))
convert(::Type{BigFloat}, x::Union{UInt8,UInt16,UInt32}) = BigFloat(convert(Culong,x))

convert(::Type{BigFloat}, x::Union{Float16,Float32}) = BigFloat(Float64(x))
convert(::Type{BigFloat}, x::Rational) = BigFloat(num(x)) / BigFloat(den(x))

function tryparse(::Type{BigFloat}, s::AbstractString, base::Int=0)
    z = BigFloat()
    err = ccall((:mpfr_set_str, :libmpfr), Int32, (Ptr{BigFloat}, Cstring, Int32, Int32), &z, s, base, ROUNDING_MODE[end])
    err == 0 ? Nullable(z) : Nullable{BigFloat}()
end

convert(::Type{Rational}, x::BigFloat) = convert(Rational{BigInt}, x)
convert(::Type{AbstractFloat}, x::BigInt) = BigFloat(x)

## BigFloat -> Integer
function unsafe_cast(::Type{Int64}, x::BigFloat, ri::Cint)
    ccall((:__gmpfr_mpfr_get_sj,:libmpfr), Cintmax_t,
          (Ptr{BigFloat}, Cint), &x, ri)
end
function unsafe_cast(::Type{UInt64}, x::BigFloat, ri::Cint)
    ccall((:__gmpfr_mpfr_get_uj,:libmpfr), Cuintmax_t,
          (Ptr{BigFloat}, Cint), &x, ri)
end

function unsafe_cast{T<:Signed}(::Type{T}, x::BigFloat, ri::Cint)
    unsafe_cast(Int64, x, ri) % T
end
function unsafe_cast{T<:Unsigned}(::Type{T}, x::BigFloat, ri::Cint)
    unsafe_cast(UInt64, x, ri) % T
end

function unsafe_cast(::Type{BigInt}, x::BigFloat, ri::Cint)
    # actually safe, just keep naming consistent
    z = BigInt()
    ccall((:mpfr_get_z, :libmpfr), Int32, (Ptr{BigInt}, Ptr{BigFloat}, Int32),
          &z, &x, ri)
    z
end
unsafe_cast(::Type{Int128}, x::BigFloat, ri::Cint) = Int128(unsafe_cast(BigInt,x,ri))
unsafe_cast(::Type{UInt128}, x::BigFloat, ri::Cint) = UInt128(unsafe_cast(BigInt,x,ri))
unsafe_cast{T<:Integer}(::Type{T}, x::BigFloat, r::RoundingMode) = unsafe_cast(T,x,to_mpfr(r))

unsafe_trunc{T<:Integer}(::Type{T}, x::BigFloat) = unsafe_cast(T,x,RoundToZero)


convert(::Type{Bool}, x::BigFloat) = x==0 ? false : x==1 ? true : throw(InexactError())
function convert(::Type{BigInt},x::BigFloat)
    isinteger(x) || throw(InexactError())
    trunc(BigInt,x)
end

function convert{T<:Integer}(::Type{T},x::BigFloat)
    isinteger(x) || throw(InexactError())
    trunc(T,x)
end

## BigFloat -> AbstractFloat
convert(::Type{Float64}, x::BigFloat) =
    ccall((:mpfr_get_d,:libmpfr), Float64, (Ptr{BigFloat},Int32), &x, ROUNDING_MODE[end])
convert(::Type{Float32}, x::BigFloat) =
    ccall((:mpfr_get_flt,:libmpfr), Float32, (Ptr{BigFloat},Int32), &x, ROUNDING_MODE[end])
# TODO: avoid double rounding
convert(::Type{Float16}, x::BigFloat) = convert(Float16, convert(Float32, x))

(::Type{Float64})(x::BigFloat, r::RoundingMode) =
    ccall((:mpfr_get_d,:libmpfr), Float64, (Ptr{BigFloat},Int32), &x, to_mpfr(r))
(::Type{Float32})(x::BigFloat, r::RoundingMode) =
    ccall((:mpfr_get_flt,:libmpfr), Float32, (Ptr{BigFloat},Int32), &x, to_mpfr(r))
# TODO: avoid double rounding
(::Type{Float16})(x::BigFloat, r::RoundingMode) =
    convert(Float16, Float32(x, r))

promote_rule{T<:Real}(::Type{BigFloat}, ::Type{T}) = BigFloat
promote_rule{T<:AbstractFloat}(::Type{BigInt},::Type{T}) = BigFloat
promote_rule{T<:AbstractFloat}(::Type{BigFloat},::Type{T}) = BigFloat
