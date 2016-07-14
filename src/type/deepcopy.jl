
function Base.deepcopy_internal(x::BigFloat, stackdict::ObjectIdDict)
    if haskey(stackdict, x)
        return stackdict[x]
    end
    N = precision(x)
    y = BigFloat(zero(Clong), zero(Cint), zero(Clong), C_NULL)
    ccall((:mpfr_init2,:libmpfr), Void, (Ptr{BigFloat}, Clong), &y, N)
    finalizer(y, cglobal((:mpfr_clear, :libmpfr)))
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Int32),
          &y, &x, ROUNDING_MODE[end])
    stackdict[x] = y
    return y
end
