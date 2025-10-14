# CalibrationTemplateFits.jl

[![Build Status](https://github.com/tdixon97/CalibrationTemplateFits/workflows/CI/badge.svg)](https://github.com/tdixon97/CalibrationTemplateFits/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/tdixon97/CalibrationTemplateFits.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tdixon97/CalibrationTemplateFits.jl)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)

Tools to model calibration data (with [BAT.jl](https://bat.github.io/BAT.jl/stable/)) to templates based on MC simulations, allowing also for morphing of templates with [Interpolations.jl](https://juliamath.github.io/Interpolations.jl/stable/) to include systematic uncertainties.

This can be used to fit the calibration data including the source position or HPGe response model as parameters of the statistical model.
