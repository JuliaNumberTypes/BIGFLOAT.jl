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

include("initialize.jl")        # module __init__ and consts

include("type/bigfloat.jl")     # precision parameterized type
include("type/precision.jl")    # get/reset var bit precision, get/set current default
include("type/eminemax.jl")     # exponent min and max
include("type/deepcopy.jl")     # <???>
include("type/cvtother.jl")

include("basics/comparison.jl") # total order relations
include("basics/fpnumerics.jl") #
include("basics/rounding.jl")   #
include("basics/floorceil.jl")  #

include("type/convert.jl")      # numeric type interconversions with BigFloat

include("math/arithmetic.jl")
include("math/power.jl")
include("math/elementary.jl")

end # module MPFR

