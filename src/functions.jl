
raw"""
    piecewise_linear_activeness(d::Real;fccd::Real,dlf::Real)->Real

Piecewise linear model which is "1" when `d` is more than `fccd`, is
0 when `d` less than `fccd*dlf` and is linear inbetween.

"""
function piecewise_linear_activeness(d::Real; fccd::Real, dlf::Real)

    if (d>fccd)
        return 1.0
    elseif (d<dlf*fccd)
        return 0
    else
        dl = fccd*dlf
        return (d-dl)/(fccd-dl)
    end
end
