using Aqua
using CalibrationTemplateFits

include("test_generalised_hist.jl")
include("test_functions.jl")
include("test_likelihood.jl")
include("test_utils.jl")
include("test_io.jl")

Aqua.test_all(CalibrationTemplateFits)

Test.@testset verbose=true "Package CalibrationTemplateFits" begin
    include("test_aqua.jl")
    include("test_generalised_hist.jl")
    include("test_functions.jl")
    include("test_likelihood.jl")
    include("test_utils.jl")
    include("test_io.jl")

end # testset
