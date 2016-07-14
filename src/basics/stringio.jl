
function string{P}(x::BigFP{P})
    # In general, the number of decimal places needed to read back the number exactly
    # is, excluding the most significant, ceil(log(10, 2^precision(x)))
    k = ceil(Int32, precision(x) * 0.3010299956639812)
    lng = k + Int32(8) # Add space for the sign, the most significand digit, the dot and the exponent
    buf = Array{UInt8}(lng + 1)
    # format strings are guaranteed to contain no NUL, so we don't use Cstring
    lng = ccall((:mpfr_snprintf,:libmpfr), Int32, (Ptr{UInt8}, Culong, Ptr{UInt8}, Ptr{BigFP{P}}...), buf, lng + 1, "%.Re", &x)
    if lng < k + 5 # print at least k decimal places
        lng = ccall((:mpfr_sprintf,:libmpfr), Int32, (Ptr{UInt8}, Ptr{UInt8}, Ptr{BigFP{P}}...), buf, "%.$(k)Re", &x)
    elseif lng > k + 8
        buf = Array{UInt8}(lng + 1)
        lng = ccall((:mpfr_snprintf,:libmpfr), Int32, (Ptr{UInt8}, Culong, Ptr{UInt8}, Ptr{BigFP{P}}...), buf, lng + 1, "%.Re", &x)
    end
    n = (1 <= x < 10 || -10 < x <= -1 || x == 0) ? lng - 4 : lng
    return String(buf[1:n])
end

print{P}(io::IO, b::BigFP{P}) = print(io, string(b))
show{P}(io::IO, b::BigFP{P}) = print(io, string(b))

print{P}(b::BigFP{P}) = print(Base.STDOUT, string(b))
show{P}(b::BigFP{P}) = print(Base.STDOUT, string(b))
