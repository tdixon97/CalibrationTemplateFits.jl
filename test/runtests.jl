using Aqua
using CalibrationTemplateFits

include("test_generalised_hist.jl")
include("test_functions.jl")
include("test_likelihood.jl")
include("test_utils.jl")
include("test_io.jl")
include("test_stats.jl")

if Sys.WORD_SIZE == 64
    Aqua.test_all(CalibrationTemplateFits)
else
    Aqua.test_all(
        CalibrationTemplateFits,
        stale_deps = (ignore = [:CairoMakie, :LegendMakie],),
        persistent_tasks = false,
    )
end

Test.@testset verbose=true "Package CalibrationTemplateFits" begin
    include("test_aqua.jl")
    include("test_generalised_hist.jl")
    include("test_functions.jl")
    include("test_likelihood.jl")
    include("test_utils.jl")
    include("test_io.jl")
    include("test_stats.jl")

end # testset
