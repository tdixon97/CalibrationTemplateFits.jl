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
export extract_grid_values

export build_likelihood
export spectrum_likelihood
export plot_posteriors
export plot_reconstruction
export plot_reconstruction_makie

export build_prior

export normalised_poisson_residual

include("generalised_histogram.jl")
include("functions.jl")
include("likelihood.jl")
include("io.jl")
include("stats.jl")

include("utils.jl")
include("plots.jl")

end # module CalibrationTemplateFits
