module CalibrationTemplateFits

export GeneralisedHistogram
export HistogramWithPars
export piecewise_linear_activeness

include("generalised_histogram.jl")
include("functions.jl")
include("likelihood.jl")
include("io.jl")

end # module CalibrationTemplateFits
