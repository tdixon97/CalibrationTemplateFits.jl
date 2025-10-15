
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
