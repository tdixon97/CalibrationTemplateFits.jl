using BAT
using Plots
using PDFmerger: append_pdf!

"""
    plot_model_predictions(path::String,samples::DensitySampleVector,data::Dict,models::Dict)
Plot the best fit model predictions
"""
function plot_model_predictions(path::String,samples::DensitySampleVector,data::Dict,models::Dict)

    
    
end

"""
    function make_summary_plots(path::String, samples::DensitySamplesVector)

Produce some summary plots of the fit results.
"""
function make_summary_plots(path::String, samples::DensitySampleVector, dets::AbstractVector, vary_fccd::Bool)

    Plots.default(
        framestyle = :box,
        background_color = :white,
        titlefontsize = 10,
        guidefontsize = 10,
        tickfontsize = 10,
        legendfontsize = 8,
    )
    
    color_scheme = [:red4, :red, :salmon]
    
    if ispath(path)
        rm(path)
    end
    
    for par in [:z, :A, :φ]
        p = plot(
            samples,
            par,
            mean = false,
            std = false,
            globalmode = false,
            marginalmode = true,
            alpha = 0.7,
            colors = color_scheme,
            nbins = 500,
            title = "Marginalized Distribution for $par",
            size = (600, 400),
        )
        
        savefig(p, "temp.pdf")
        append_pdf!(path, "temp.pdf", cleanup = true)
    end

    if (vary_fccd)

        for det in dets
            par = Symbol("fccd_$det")
            p = plot(
            samples,
            Symbol("fccd_$det"),
            mean = false,
            std = false,
            globalmode = false,
            marginalmode = true,
            alpha = 0.7,
            colors = color_scheme,
            nbins = 500,
            title = "Marginalized Distribution for $par",
            size = (600, 400),
            )
        
            savefig(p, "temp.pdf")
            append_pdf!(path, "temp.pdf", cleanup = true)
        end
    end

    
end
