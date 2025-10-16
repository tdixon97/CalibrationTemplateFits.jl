using BAT
using Plots
using PDFmerger: append_pdf!
default(
    framestyle = :box,               # Grid line transparency
    background_color = :white,       # Background color of the plot,
    titlefontsize = 10,     # Global title font size
    guidefontsize = 10,     # Global axis label font size
    tickfontsize = 10,      # Global tick label font size
    legendfontsize = 8,     # Global legend font size
)


tol_colors = ColorSchemes.tol_muted
color_schemes = [:red4, :red, :salmon]

"""
    function make_summary_plots(path::String, samples::DensitySamplesVector)

Produce some summary plots of the fit results.
"""
function make_summary_plots(path::String, samples::DensitySamplesVector)


    color_scheme = [:red4, :red, :salmon]

    for par in [:z,:A,:φ]
        p = plot(
            samples, par,
            mean = false, std = false, globalmode = false, marginalmode = true, alpha=0.7,
                                colors = color_scheme,
            nbins = 200, title = "Marginalized Distribution for z offset", size=(800,500)
        )

        savefig(p, "temp.pdf")
        append_pdf(path,
                    "temp.pdf",
                    cleanup = true,
                )
    end
end
end
