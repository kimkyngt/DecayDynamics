using Distributions, Roots

# Define a custom distribution type
struct DoubleExp{T<:Real} <: ContinuousUnivariateDistribution
    tauP::T
    tauD::T
end

# Implement the methods required by the ContinuousUnivariateDistribution interface

# Probability density function (PDF)
function Distributions.pdf(d::DoubleExp, x::Real)
    if x >= 0
        return 1/(d.tauP - d.tauD)*(exp(-x/d.tauP) - exp(-x/d.tauD))
    else
        return 0
    end
end

# Cumulative distribution function (CDF)
function Distributions.cdf(d::DoubleExp, x::Real)
    if x >= 0
        return 1 + 1/(d.tauP - d.tauD)*(-d.tauP*exp(-x/d.tauP) + d.tauD*exp(-x/d.tauD))
    else
        return 0
    end
end

# Quantile function (inverse CDF)
function Distributions.quantile(d::DoubleExp, p::Real)
    @assert 0 ≤ p ≤ 1 "p must be between 0 and 1"
    f(x) = 1 + 1/(d.tauP - d.tauD)*(-d.tauP*exp(-x/d.tauP) + d.tauD*exp(-x/d.tauD)) - p
    return find_zero(f, (0, Inf))
end

# Support of the distribution
function Distributions.support(d::DoubleExp)
    return Interval(0, Inf)
end

# Mean of the distribution
function Distributions.mean(d::DoubleExp)
    return d.tauP*d.tauD/(d.tauD - d.tauP)
end

# Variance of the distribution
function Distributions.var(d::DoubleExp)
    return d.tauP^2 * d.tauD^2 / ((d.tauD - d.tauP)^2 * (d.tauD + d.tauP))
end