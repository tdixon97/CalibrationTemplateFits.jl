module CalibrationTemplateFits

export GeneralisedHistogram
export HistogramWithPars
export piecewise_linear_activeness
export read_data_histograms
export make_data_histograms

export parse_binning
export read_models
export read_models_evt
export read_models_hist

export build_likelihood
export spectrum_likelihood
export make_summary_plots
export build_prior

include("generalised_histogram.jl")
include("functions.jl")
include("likelihood.jl")
include("io.jl")
include("utils.jl")
include("plots.jl")

end # module CalibrationTemplateFits
