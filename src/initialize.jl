
function __init__()
    try
        # set exponent to full range by default
        set_emin!(get_emin_min())
        set_emax!(get_emax_max())
    catch ex
        Base.showerror_nostdio(ex,
            "WARNING: Error during initialization of module MPFR")
    end
end

const ROUNDING_MODE = Cint[0]
const BIGFLOAT_PRECISION = [256]
