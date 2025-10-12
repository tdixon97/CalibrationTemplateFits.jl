using Aqua
using CalibrationTemplateFits

include("test_generalised_hist.jl")

Aqua.test_all(CalibrationTemplateFits)

Test.@testset verbose=true "Package CalibrationTemplateFits" begin
    include("test_aqua.jl")
    include("test_generalised_hist.jl")
end # testset
