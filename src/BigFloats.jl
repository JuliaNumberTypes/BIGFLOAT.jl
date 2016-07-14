# This file is a part of Julia. License is MIT: http://julialang.org/license

module MPFR

export
    BigFloat,
    setprecision,
    big_str

import
    Base: (*), +, -, /, <, <=, ==, >, >=, ^, besselj, besselj0, besselj1, bessely,
        bessely0, bessely1, ceil, cmp, convert, copysign, div,
        exp, exp2, exponent, factorial, floor, fma, hypot, isinteger,
        isfinite, isinf, isnan, ldexp, log, log2, log10, max, min, mod, modf,
        nextfloat, prevfloat, promote_rule, rem, round, show,
        sum, sqrt, string, print, trunc, precision, exp10, expm1,
        gamma, lgamma, digamma, erf, erfc, zeta, eta, log1p, airyai,
        eps, signbit, sin, cos, tan, sec, csc, cot, acos, asin, atan,
        cosh, sinh, tanh, sech, csch, coth, acosh, asinh, atanh, atan2,
        cbrt, typemax, typemin, unsafe_trunc, realmin, realmax, rounding,
        setrounding, maxintfloat, widen, significand, frexp, tryparse

import Base.Rounding: rounding_raw, setrounding_raw
import Base.GMP: ClongMax, CulongMax, CdoubleMax, Limb
import Base.Math.lgamma_r

include("type/type_mpfr.jl")



end # module MPFR

