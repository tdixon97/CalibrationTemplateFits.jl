using BAT
using Plots
using PDFmerger: append_pdf!

Plots.default(
    framestyle = :box,
    background_color = :white,
    titlefontsize = 10,
    guidefontsize = 10,
    tickfontsize = 10,
    legendfontsize = 8,
)

color_schemes = [:red4, :red, :salmon]

"""
    function make_summary_plots(path::String, samples::DensitySamplesVector)

Produce some summary plots of the fit results.
"""
function make_summary_plots(path::String, samples::DensitySampleVector)

    color_scheme = [:red4, :red, :salmon]

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
            nbins = 200,
            title = "Marginalized Distribution for $par",
            size = (600, 400),
        )

        savefig(p, "temp.pdf")
        append_pdf!(path, "temp.pdf", cleanup = true)
    end
end
