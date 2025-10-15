
using TypedTables
using CalibrationTemplateFits
@testset "test_read_data" begin

    path = joinpath(@__DIR__, "test_files")

    @test CalibrationTemplateFits.read_data(path, ["evt_data.lh5"]) isa Table

end
