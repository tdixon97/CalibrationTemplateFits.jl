
using TypedTables
using CalibrationTemplateFits
using StatsBase

@testset "test_read_data" begin

    path = joinpath(@__DIR__, "test_files")

    @test CalibrationTemplateFits.read_data(path, ["evt_data.lh5"]) isa Table

    hists = make_data_histograms(
        path,
        ["evt_data.lh5"],
        Dict(:det1=>1, :det2=>2, :det3=>3),
        0:10:4000,
    )
    @test all(x isa Histogram for x in values(hists))
    @test Set(keys(hists)) == Set([:det1, :det2, :det3])

    # should also work with vector
    hists = make_data_histograms(
        path,
        ["evt_data.lh5"],
        Dict(:det1=>1, :det2=>2, :det3=>3),
        [2500, 2600, 2620, 2700],
    )
    @test all(x isa Histogram for x in values(hists))
    @test Set(keys(hists)) == Set([:det1, :det2, :det3])

end

@testset "test_read_data_hist" begin

    path = joinpath(@__DIR__, "test_files")

    hists = read_data_histograms(path*"pdf_data.lh5", [:det1, :det2, :det3], 0:10:4000)
    @test all(x isa Histogram for x in values(hists))
    @test Set(keys(hists)) == Set([:det1, :det2, :det3])

end

@testset "test_read_mc" begin
    path = joinpath(@__DIR__, "test_files")


    @test CalibrationTemplateFits.read_mc_files("det1", [String(path)*"/mc.lh5"]) isa Vector
    @test CalibrationTemplateFits.extract_mc_coords(
        "z_-1_phi_0",
        r".*z_([-\d.]+)_phi_([-\d.]+)",
    ) == (-1, 0)


    models = read_models(
        [:det1, :det2, :det3],
        [
            path*"/hit_files/z_-1_phi_0",
            path*"/hit_files/z_1_phi_0",
            path*"/hit_files/z_-1_phi_1",
            path*"/hit_files/z_1_phi_1",
        ],
        2600:20:2620,
        r".*z_([-\d.]+)_phi_([-\d.]+)",
    )

    @test all(mod isa GeneralisedHistogram for mod in values(models))
end

@testset "test_read_histograms" begin

    path = joinpath(@__DIR__, "test_files")
    files = [
        path*"/pdf_files/z_-1_phi_0.lh5",
        path*"/pdf_files/z_1_phi_0.lh5",
        path*"/pdf_files/z_-1_phi_1.lh5",
        path*"/pdf_files/z_1_phi_1.lh5",
    ]

    @test CalibrationTemplateFits.read_hist("hit/det1", files[1]) isa Histogram
    mods = read_models_hist(
        ["det1", "det2", "det3"],
        files,
        2600:30:2630,
        r".*z_([-\d.]+)_phi_([-\d.]+)",
        "hit",
    )
    @test all(mod isa GeneralisedHistogram for mod in values(mods))
end
