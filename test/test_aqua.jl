import Test
import Aqua
import CalibrationTemplateFits

Test.@testset "Package ambiguities" begin
    Test.@test isempty(Test.detect_ambiguities(CalibrationTemplateFits))
end # testset

Test.@testset "Aqua tests" begin
    Aqua.test_all(CalibrationTemplateFits, ambiguities = true)
end # testset
