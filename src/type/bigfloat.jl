# Basic type and initialization definitions

type BigFP{P} <: AbstractFloat
    prec::Clong
    sign::Cint
    exp::Clong
    d::Ptr{Limb}
    
    function BigFP{P}()
        N = precision(BigFP)
        z = new(zero(Clong), zero(Cint), zero(Clong), C_NULL)
        ccall((:mpfr_init2,:libmpfr), Void, (Ptr{BigFP}, Clong), &z, N)
        finalizer(z, cglobal((:mpfr_clear, :libmpfr)))
        return z
    end
    # Not recommended for general use
    function BigFP{P}(prec::Clong, sign::Cint, exp::Clong, d::Ptr{Void})
        new(prec, sign, exp, d)
    end
end
