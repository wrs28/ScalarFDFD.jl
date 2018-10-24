################################################################################
### PML
################################################################################
const EXTINCTION = 1e-8
const POWER_LAW = 2
const SCALING_ANGLE = .25


################################################################################

const MINIMUM_F = 1e-16

const TRANSVERSE_MODE_RESOLUTION_MULTIPLIER = 3

const PROPAGATION_CONSTANT_IMAGINARY_TOLERANCE = 1e-2

const NUMBER_OF_QUADRATURE_POINTS_ANG_MOM_ANALYSIS = 1001

const PROGRESS_UPDATE_TIME = .05 # in seconds

const DEFAULT_SUBSAMPLE_NUMBER = 20 # default number of samples used in sub-pixel smoothing (total number is square of this)


########### FOR PLOTTING
function __init__()

    if !haskey(ENV, "SCALAR_FDFD_COLOR_THEME")
        SCALAR_FDFD_COLOR_THEME = :default
    else
        SCALAR_FDFD_COLOR_THEME = Symbol(ENV["SCALAR_FDFD_COLOR_THEME"])
    end

    global DEFAULT_COLOR_SCHEME = SCALAR_FDFD_COLOR_THEME
end
