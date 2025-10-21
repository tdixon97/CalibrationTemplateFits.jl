using BAT
using Plots
using StatsBase
using PDFmerger: append_pdf!
using Unitful
using Measures
"""
    function plot_reconstruction(data::Dict,models::Dicts,out_path::String, mode::NamedTuple)

Produce some summary plots of the fit results.
"""
function plot_reconstruction(
    data::Dict,
    models::Dict,
    out_path::String,
    mode::NamedTuple,
    norm::Float64,
)

    Plots.default(
        framestyle = :box,
        background_color = :white,
        titlefontsize = 10,
        guidefontsize = 10,
        tickfontsize = 10,
        legendfontsize = 8,
    )

    # number of dets
    n = length(keys(data))
    dets = []

    hData = fit(Histogram{Float64}, Float64[], 0.5:1:(length(keys(data))+0.5))
    hPred = fit(Histogram{Float64}, Float64[], 0.5:1:(length(keys(data))+0.5))

    # loop
    for (idx, det) in enumerate(keys(data))

        hData.weights[idx] = sum(data[det].weights)
        hPred.weights[idx] = sum(get_weights(models[det]; mode...))*mode.A*norm
        push!(dets, det)
    end

    # plot the data
    p = plot(
        hData,
        st = :steps,
        yscale = :log10,
        ylabel = "Counts",
        label = "Data",
        linewidth = 2,
        size = (1200, 400),
        tickfont = font(10),
        xrotation = 90,
        bottommargin = 3mm,
        leftmargin = 15mm,
        ylims = (0.01, maximum(hData.weights)*1.5),
    )

    # plot the best fit
    plot!(
        hPred,
        st = :steps,
        yscale = :log10,
        ylabel = "Counts",
        label = "Best fit",
        linewidth = 2,
        size = (1400, 400),
        tickfont = font(10),
        xrotation = 90,
        bottommargin = 3mm,
        leftmargin = 15mm,
        ylims = (0.01, maximum(hData.weights)*1.5),
        xticks = false,
    )

    # add residual
    xlims!(0.01, n + 1 + 0.01)
    y1, y2 = ylims(p)

    r = normalised_poisson_residual.(hPred.weights, hData.weights)

    p2 = scatter(
        1:n,
        r,
        rotation = 90,
        marker = ".",
        label = nothing,
        size = (800, 400),
        ylims = (-6, 6),
        xticks = (1:n, dets),
        tickfont = font(10),
        bottommargin = 15mm,
        leftmargin = 10mm,
        color = :black,
        xlims = (0.01, n+0.5),
    )

    # and some pretty bands
    x1, x2 = xlims(p2)
    y1, y2 = -3, 3
    plot!(
        p2,
        [x1, x2, x2, x1],
        [y1, y1, y2, y2],
        color = "red",
        alpha = 0.2,
        seriestype = :shape,
        linewidth = 0,
        label = nothing,
    )

    y1, y2 = -2, 2
    plot!(
        p2,
        [x1, x2, x2, x1],
        [y1, y1, y2, y2],
        color = "orange",
        alpha = 0.4,
        seriestype = :shape,
        linewidth = 0,
        label = nothing,
    )

    y1, y2 = -1, 1
    plot!(
        p2,
        [x1, x2, x2, x1],
        [y1, y1, y2, y2],
        color = "green",
        alpha = 0.3,
        seriestype = :shape,
        linewidth = 0,
        label = nothing,
    )

    scatter!(
        1:n,
        r,
        rotation = 90,
        marker = ".",
        label = nothing,
        size = (1000, 400),
        ylims = (-maximum(abs.(r))-1, maximum(abs.(r))+1),
        xticks = (1:n, dets),
        tickfont = font(10),
        bottommargin = 15mm,
        leftmargin = 10mm,
        color = :black,
    )
    ylabel!("Residual [σ]")

    p3=plot(
        p,
        p2,
        framestyle = :box,
        layout = Plots.grid(2, 1, heights = [0.70, 0.3]),
        size = (1000, 500),
        link = :x,
        bottom_margin = 9mm,
        left_margin = 5mm,
        right_margin = 5mm,
    )

    # save
    savefig(p3, out_path)



end

"""
    function plot_posteriors(path::String, samples::DensitySamplesVector)

Produce some summary plots of the fit results.
"""
function plot_posteriors(
    path::String,
    samples::DensitySampleVector,
    dets::AbstractVector,
    vary_fccd::Bool,
)

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
