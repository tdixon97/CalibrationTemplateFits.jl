using BAT
using Plots
using StatsBase
using PDFmerger: append_pdf!
using Unitful
using Measures
using CairoMakie
using LegendMakie
using StatsBase
using PDFmerger: append_pdf!
using FilePathsBase: ispath
using CairoMakie
using StatsBase
using PDFmerger: append_pdf!

"""
Simple way to plot a histogram as bar.
"""
function plot_hist!(ax, h::Histogram; kwargs...)

    counts, bins = h.weights, h.edges[1]
    counts_eps = counts .+ 1e-5
    bins = collect(bins)

    append!(bins, bins[end]+diff(bins)[end])
    append!(counts_eps, [1e-5])

    barplot!(
        ax,
        bins[1:(end-1)]+diff(bins)/2,
        counts_eps,
        gap = 0.0,
        width = diff(bins);
        kwargs...,
    )

end

function plot_reconstruction_makie(
    data::Dict,
    models::Dict,
    out_path::String,
    mode::NamedTuple,
    norm::Float64,
)
    if ispath(out_path)
        rm(out_path)
    end

    dets = collect(keys(data))
    n = length(dets)

    # Create aggregate histograms
    hData = fit(Histogram{Float64}, Float64[], 0.5:1:(n+0.5))
    hPred = fit(Histogram{Float64}, Float64[], 0.5:1:(n+0.5))

    for (idx, det) in enumerate(dets)
        hData.weights[idx] = sum(data[det].weights)
        hPred.weights[idx] = sum(get_weights(models[det]; mode...)) * mode.A * norm
    end

    # --- Global summary plot ---
    fig = Figure(size = (1000, 500))

    ax1 = Axis(
        fig[1, 1];
        ylabel = "Counts",
        yscale = log10,
        title = "Global fit summary",
        xticklabelsvisible = false,
    )

    # Use hist! to draw bar-style histograms
    plot_hist!(ax1, hData, label = "Data", alpha = 0.4)
    w = append!([1e-5], hPred.weights .+= 1e-5)
    append!(w, [1e-5])

    e = collect(hPred.edges[1])
    append!(e, e[end]+diff(e)[end])

    CairoMakie.stairs!(ax1, e, w, color = "orange", label = "Best fit")

    CairoMakie.ylims!(ax1, 0.01, maximum(hData.weights) * 1.5)
    axislegend(ax1; position = :rt)

    # Residuals subplot
    ax2 = Axis(fig[2, 1]; ylabel = "Residual [σ]", xticks = (1:n, dets))

    r = normalised_poisson_residual.(hPred.weights, hData.weights)
    band!(ax2, 0:n, -3*ones(n+1), 3*ones(n+1), color = (:turquoise2, 0.3))
    lines!(ax2, 0:n, zeros(n+1), color = :blue, linewidth = 2)

    CairoMakie.scatter!(ax2, 1:n, r, color = :black, markersize = 10)

    ax2.xticklabelrotation[] = 90
    ax2.xticklabelsize[] = 12

    linkxaxes!(ax1, ax2)

    save("temp.pdf", fig)
    append_pdf!(out_path, "temp.pdf", cleanup = true)

    # --- Individual detector plots ---
    if length(data[dets[1]].weights) == 1
        return
    end

    for det in dets
        hData = data[det]
        hPred = get_histogram(models[det]; mode...)
        hPred.weights .*= mode.A * norm

        xl = (first(hData.edges[1]), last(hData.edges[1]))
        r = normalised_poisson_residual.(hPred.weights, hData.weights)

        fig = Figure(size = (1000, 500))

        ax1 = Axis(
            fig[1, 1];
            ylabel = "Counts",
            title = det,
            yscale = log10,
            xticklabelsvisible = false,
        )

        plot_hist!(ax1, hData, label = "Data", alpha = 0.4)
        w = append!([1e-5], hPred.weights .+= 1e-5)
        append!(w, [1e-5])

        e = collect(hPred.edges[1])
        append!(e, e[end]+diff(e)[end])

        CairoMakie.stairs!(ax1, e, w, color = "orange", label = "Best fit")

        CairoMakie.ylims!(ax1, 0.01, maximum(hData.weights) * 1.5)
        CairoMakie.xlims!(ax1, xl...)
        axislegend(ax1; position = :rt)

        ax2 = Axis(fig[2, 1]; ylabel = "Residual [σ]", xlabel = "Energy [keV]")

        CairoMakie.ylims!(ax2, -6, 6)

        band!(ax2, e, -3*ones(length(e)), 3*ones(length(e)), color = (:darkcyan, 0.3))
        lines!(ax2, e, zeros(length(e)), color = :blue, linewidth = 2)

        CairoMakie.scatter!(
            ax2,
            hData.edges[1][1:(end-1)],
            r,
            color = :black,
            markersize = 10,
        )

        linkxaxes!(ax1, ax2)

        save("temp.pdf", fig)
        append_pdf!(out_path, "temp.pdf", cleanup = true)
    end
end


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
    if ispath(out_path)
        rm(out_path)
    end
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
    xlims!(0.5, n + 0.5)
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
        xlims = (0.5, n + 0.5),
    )

    # and some pretty bands
    x1, x2 = xlims(p2)
    y1, y2 = -3, 3
    plot!(
        p2,
        [x1, x2, x2, x1],
        [y1, y1, y2, y2],
        color = "turquoise2",
        alpha = 0.4,
        seriestype = :shape,
        linewidth = 0,
        label = "3σ",
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
    savefig(p3, "temp.pdf")
    append_pdf!(out_path, "temp.pdf", cleanup = true)

    # now plot the overall reconstruction - only if we have multiple bins
    if (length(data[dets[1]].weights)==1)
        return
    end

    # loop over detectors
    for det in dets


        hData = data[det]
        hPred = get_histogram(models[det]; mode...)
        hPred.weights*=mode.A*norm

        xl = (collect(hData.edges)[1], collect(hData.edges)[end])
        # plot the data
        p = plot(
            hData,
            st = :steps,
            ylabel = "Counts",
            label = "Data",
            linewidth = 2,
            size = (1200, 400),
            tickfont = font(10),
            xrotation = 90,
            bottommargin = 3mm,
            title = det,
            leftmargin = 15mm,
            ylims = (0.01, maximum(hData.weights)*1.2),
            xlims = xl,
        )

        # plot the best fit
        plot!(
            hPred,
            st = :steps,
            ylabel = "Counts",
            label = "Best fit",
            linewidth = 2,
            size = (1400, 400),
            tickfont = font(10),
            xrotation = 90,
            bottommargin = 3mm,
            leftmargin = 15mm,
            ylims = (0.01, maximum(hData.weights)*1.2),
            xlims = xl,
        )

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
            xlims = xl,
        )

        # and some pretty bands
        x1, x2 = xlims(p2)
        y1, y2 = -3, 3

        plot!(
            p2,
            [x1, x2, x2, x1],
            [y1, y1, y2, y2],
            color = "turquoise2",
            alpha = 0.4,
            seriestype = :shape,
            linewidth = 0,
            label = "3σ",
        )

        scatter!(
            1:n,
            r,
            rotation = 90,
            marker = ".",
            label = nothing,
            size = (1000, 400),
            ylims = (-maximum(abs.(r))-1, maximum(abs.(r))+1),
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
        savefig(p3, "temp.pdf")
        append_pdf!(out_path, "temp.pdf", cleanup = true)


    end


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
        p = Plots.plot(
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
            p = Plots.plot(
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
