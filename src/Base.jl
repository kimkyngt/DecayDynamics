# Base package for the project
using QuantumOptics, WignerSymbols, LinearAlgebra

""" 
Get |F, m⟩
"""
function Fm_state(F, m)
    if abs(m) > F
        throw(ArgumentError("m must be in the range [-F, F]"))
    end
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
function σ(F_1::Number, m_1::Number, F_2::Number, m_2::Number)
    ψe = Fm_state(F_1, m_1)⊕ Ket(SpinBasis(F_2))
    ψg = Ket(SpinBasis(F_1)) ⊕ Fm_state(F_2, m_2)
    return sparse(ψe ⊗ dagger(ψg))
end


"""
Σ_q(q::Int, F_i::Array, k::Int=1, l::Int=2)
Atomic lowering operator with hyperfine structure. Method for multiple levels.
- `q`: -1, 0, 1
- `F_i`: Vector of the spin of each level.
- `k`: index for the upper state.
- `l`: index for the lower state.
"""
function Σ_q(q::Int, F_i::Array, k::Int=1, l::Int=2)
    if abs(q) > 1
        throw(ArgumentError("Argument q must be one of [-1, 0, 1]"))
    end
    
    embed(
        directsum([SpinBasis(F) for F in F_i]...),
        directsum([SpinBasis(F) for F in F_i]...),
        [k, l],
        dagger(sum(clebschgordan(F_i[l], m, 1, q, F_i[k])*σ(F_i[k], m + q, F_i[l], m) for m = -F_i[l]:1:F_i[l]if abs(m + q) <= F_i[k] ))
    )
end

"""
Σ_iq(i::Int, q::Int, F_i::Vector{Rational} ; kindx::Int=1).
Hyperfine lowering operator for multilevel two atoms. 
- `i`: index of the atom under interest.
- `q`: -1, 0, 1
- `F_i`: Array of the spin of each level.
- `k`: index for the upper state.
- `l`: index for the lower state.
"""
function Σ_iq(i::Int, q::Int, F_i::Vector, k::Int, l::Int)
    if i > 2
        throw(ArgumentError("Argument i must be one of [1, 2]"))
    end
    if abs(q) > 1
        throw(ArgumentError("Argument q must be one of [-1, 0, 1]"))
    end
    sigma_single = Σ_q(q, F_i, k, l)
    embed(sigma_single.basis_l ⊗ sigma_single.basis_l, [i], [sigma_single])
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
Modified from CollectiveSpins.jl package. 
Static limit, imaginary part goes to 2/3
"""
function GreenTensor(r::Vector{<:Number},k::Real=2π)
    n = norm(r)
    if n == 0
        return sparse(im*Matrix(I,3,3))*2/3
    else
        rn = r./n
        return exp(im*k*n)*(
            (1/(k*n) + im/(k*n)^2 - 1/(k*n)^3) .* Matrix(I,3,3) +
            -(1/(k*n) + 3im/(k*n)^2 - 3/(k*n)^3) .* (rn*rn')
        )
    end
end