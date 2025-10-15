
using TypedTables
using CalibrationTemplateFits
using StatsBase

@testset "test_read_data" begin

    path = joinpath(@__DIR__, "test_files")

    @test CalibrationTemplateFits.read_data(path, ["evt_data.lh5"]) isa Table

    hists = read_data_histograms(
        path,
        ["evt_data.lh5"],
        Dict(:det1=>1, :det2=>2, :det3=>3),
        0:10:4000,
    )
    @test all(x isa Histogram for x in values(hists))
    @test Set(keys(hists)) == Set([:det1, :det2, :det3])

    # should also work with vector
    hists = read_data_histograms(
        path,
        ["evt_data.lh5"],
        Dict(:det1=>1, :det2=>2, :det3=>3),
        [2500, 2600, 2620, 2700],
    )
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
        [path*"/z_-1_phi_0", path*"/z_1_phi_0", path*"/z_-1_phi_1", path*"/z_1_phi_1"],
        2600:20:2620,
        r".*z_([-\d.]+)_phi_([-\d.]+)",
    )

    @test all(mod isa GeneralisedHistogram for mod in values(models))
end
