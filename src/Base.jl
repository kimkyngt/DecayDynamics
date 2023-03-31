# Base package for the project
using QuantumOptics, WignerSymbols, LinearAlgebra

""" 
Get |F, m>
"""
function Fm_state(F, m)
    b = SpinBasis(F)
    return normalize(sigmam(b)^Int(F - m) *spinup(b))
end

"""
spherical basis vectors
"""
function ϵ_q(q::Int)
    if q == 0
        return transpose([0, 0, 1])
    elseif abs(q) == 1
        return -q/sqrt(2) * transpose([1, q*im, 0])
    else 
        throw(ArgumentError("Argument q must be one of [-1, 0, 1]"))
    end
end

""" 
Atomic coherence operator.
"""
function σ(F_1, m_1, F_2, m_2)
    ψe = Fm_state(F_1, m_1)⊕ Ket(SpinBasis(F_2))
    ψg = Ket(SpinBasis(F_1)) ⊕ Fm_state(F_2, m_2)
    return sparse(ψe ⊗ dagger(ψg))
end


"""
Atomic lowering operator with hyperfine structure. Spherical basis.
"""
function Σ_q(q::Int, F_1, F_2)
    return sparse(dagger(sum(clebschgordan(F_2, m, 1, q, F_1)*σ(F_1, m + q, F_2, m) for m = -F_2:1:F_2 if abs(m + q) <= F_1 )))
end


"""
Convert polar coordinate to Cartesian for plotting.
"""
function polar_2_xyz(θ, ϕ, r=1)
    x = [sin(θi)*cos(ϕj) for θi ∈ θ,  ϕj ∈ ϕ].*r
    y = [sin(θi)*sin(ϕj) for θi ∈ θ,  ϕj ∈ ϕ].*r
    z = [cos(θi) 	     for θi ∈ θ,  ϕj ∈ ϕ].*r
    return x, y, z
end

"""
Draw a unit sphere. Nsample determines the density
"""
function draw_unit_sphere(Nsample=20, )
    θ = range(0, π, length = Nsample)
    ϕ = range(0, 2π, length  = Nsample)
    fig = plot(
        polar_2_xyz(θ, ϕ, ones(length(θ))),
        st=:surface,
        color=:grey, 
        legend=false,
        )
    return fig
end


"""
    GreenTensor(r::Vector, k::Number=2π)
Calculate the Green's Tensor at position r for wave number k defined by
```math
G = e^{ikr}\\Big[\\left(\\frac{1}{kr} + \\frac{i}{(kr)^2} - \\frac{1}{(kr)^3}\\right)*I -
    \\textbf{r}\\textbf{r}^T\\left(\\frac{1}{kr} + \\frac{3i}{(kr)^2} - \\frac{3}{(kr)^3}\\right)\\Big]
```
Choosing `k=2π` corresponds to the position `r` being given in units of the
wavelength associated with the dipole transition.
Returns a 3×3 complex Matrix.

Modified from CollectiveSpins.jl package
"""
function GreenTensor(r::Vector{<:Number},k::Real=2π)
    n = norm(r)
    if n == 0
        return sparse(im*Matrix(I,3,3))
    else
        rn = r./n
        return exp(im*k*n)*(
            (1/(k*n) + im/(k*n)^2 - 1/(k*n)^3) .* Matrix(I,3,3) +
            -(1/(k*n) + 3im/(k*n)^2 - 3/(k*n)^3) .* (rn*rn')
        )
    end
end