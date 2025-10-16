module CalibrationTemplateFits

export GeneralisedHistogram
export HistogramWithPars
export piecewise_linear_activeness
export read_data_histograms
export parse_binning
export read_models
export build_likelihood
export spectrum_likelihood
export make_summary_plots

include("generalised_histogram.jl")
include("functions.jl")
include("likelihood.jl")
include("io.jl")
include("utils.jl")
include("plots.jl")

end # module CalibrationTemplateFits
