
function Base.deepcopy_internal(x::BigFP, stackdict::ObjectIdDict)
    if haskey(stackdict, x)
        return stackdict[x]
    end
    N = precision(x)
    y = BigFP(zero(Clong), zero(Cint), zero(Clong), C_NULL)
    ccall((:mpfr_init2,:libmpfr), Void, (Ptr{BigFP}, Clong), &y, N)
    finalizer(y, cglobal((:mpfr_clear, :libmpfr)))
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{BigFP}, Ptr{BigFP}, Int32),
          &y, &x, ROUNDING_MODE[end])
    stackdict[x] = y
    return y
end
