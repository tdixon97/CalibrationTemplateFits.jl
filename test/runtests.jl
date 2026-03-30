using Aqua
using CalibrationTemplateFits

include("test_generalised_hist.jl")
include("test_functions.jl")
include("test_likelihood.jl")
include("test_utils.jl")
include("test_io.jl")
include("test_stats.jl")

Aqua.test_all(
    CalibrationTemplateFits,
    stale_deps = (ignore = [:CairoMakie, :LegendMakie],),
    persistent_tasks = false,
)

Test.@testset verbose=true "Package CalibrationTemplateFits" begin
    include("test_aqua.jl")
    include("test_generalised_hist.jl")
    include("test_functions.jl")
    include("test_likelihood.jl")
    include("test_utils.jl")
    include("test_io.jl")
    include("test_stats.jl")
    include("test_plots.jl")

end # testset
