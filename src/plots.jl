using StatsBase
using PDFmerger: append_pdf!
using FilePathsBase: ispath
using CairoMakie
using LegendMakie

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
