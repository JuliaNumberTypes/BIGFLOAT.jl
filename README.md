# MPFR3.jl
## A precision parameterized rewrite of base/mpfr.jl as the type `BigFloat`

####<p align="right">`Jeffrey Sarnoff © 2016˗Aug˗10 ≏ New York City`</p>

__there are warnings when loading this because of overwriting defs__


The tests run (or most do), some div tests with big"1" inputs are commented out because the big"1" stuff is not defined within mpfr and I don't know what is going on with them that fails.


I left mixed precision ops promoting to the smaller precision just to match behavior .. imo that is the less helpful, albeit more vanilla, choice.

#### Use


```julia
Pkg.clone("https://github.com/JuliaNumberTypes/MPFR3.jl")
using MPFR3
BigFloat = MPFR3.BigFloat

a = atan(BigFloat(1)) * 4

b = atan(BigFloat{53}(1)) * 4

c = BigFloat{53}(a)

d = (a + b) / 2

```


