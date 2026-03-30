using CalibrationTemplateFits
using Test
using StatsBase

@static if Sys.WORD_SIZE == 64
    import CairoMakie
end

@testset "test_plot_hist" begin
    @static if Sys.WORD_SIZE == 64
        fig = CairoMakie.Figure()
        ax = CairoMakie.Axis(fig[1, 1])
        h = fit(Histogram{Float64}, [1.0, 2.0, 3.0], 0.0:1.0:4.0)
        @test_nowarn plot_hist!(ax, h)
    end
end

@testset "test_plot_reconstruction" begin
    h1 = append!(Histogram(2600:100:2700), [2610])
    h2 = append!(Histogram(2600:100:2700), [2690])
    hists = [HistogramWithPars(h1, par = 0), HistogramWithPars(h2, par = 1)]
    model = GeneralisedHistogram(hists, par = 0:1:1)

    data = fit(Histogram{Float64}, [2650.0], 2600:100:2700)
    data_dict = Dict(:det1 => data)
    models_dict = Dict(:det1 => model)
    mode = (par = 0.0, A = 1.0)

    mktempdir() do tmpdir
        out_path = joinpath(tmpdir, "test_output.pdf")
        cd(tmpdir) do
            @test_nowarn plot_reconstruction(data_dict, models_dict, out_path, mode, 1.0)
        end
        @test isfile(out_path)
    end
end

@static if Sys.WORD_SIZE == 64
    @testset "test_plot_reconstruction_makie" begin
        h1 = append!(Histogram(2600:100:2700), [2610])
        h2 = append!(Histogram(2600:100:2700), [2690])
        hists = [HistogramWithPars(h1, par = 0), HistogramWithPars(h2, par = 1)]
        model = GeneralisedHistogram(hists, par = 0:1:1)

        data = fit(Histogram{Float64}, [2650.0], 2600:100:2700)
        data_dict = Dict(:det1 => data)
        models_dict = Dict(:det1 => model)
        mode = (par = 0.0, A = 1.0)

        mktempdir() do tmpdir
            out_path = joinpath(tmpdir, "test_output_makie.pdf")
            cd(tmpdir) do
                @test_nowarn plot_reconstruction_makie(
                    data_dict,
                    models_dict,
                    out_path,
                    mode,
                    1.0,
                )
            end
            @test isfile(out_path)
        end
    end
end
