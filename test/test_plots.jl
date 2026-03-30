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

# Shared test data for reconstruction tests
_h1 = append!(Histogram(2600:100:2700), [2610])
_h2 = append!(Histogram(2600:100:2700), [2690])
_hists = [HistogramWithPars(_h1, par = 0), HistogramWithPars(_h2, par = 1)]
_model = GeneralisedHistogram(_hists, par = 0:1:1)
_data = fit(Histogram{Float64}, [2650.0], 2600:100:2700)
_data_dict = Dict(:det1 => _data)
_models_dict = Dict(:det1 => _model)
_mode = (par = 0.0, A = 1.0)

@testset "test_plot_reconstruction" begin
    mktempdir() do tmpdir
        out_path = joinpath(tmpdir, "test_output.pdf")
        cd(tmpdir) do
            @test_nowarn plot_reconstruction(_data_dict, _models_dict, out_path, _mode, 1.0)
        end
        @test isfile(out_path)
    end
end

@static if Sys.WORD_SIZE == 64
    @testset "test_plot_reconstruction_makie" begin
        mktempdir() do tmpdir
            out_path = joinpath(tmpdir, "test_output_makie.pdf")
            cd(tmpdir) do
                @test_nowarn plot_reconstruction_makie(
                    _data_dict,
                    _models_dict,
                    out_path,
                    _mode,
                    1.0,
                )
            end
            @test isfile(out_path)
        end
    end
end
