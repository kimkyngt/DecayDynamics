using Plots, QuantumOptics, CollectiveSpins, LaTeXStrings, WignerSymbols

""" 
Get |F, m>
"""
function Fm_state(F::Rational, m::Rational)
    b = SpinBasis(F)
    return normalize(sigmam(b)^Int(F - m) *spinup(b))
end

""" 
Atomic coherence operator for e-g system.
"""
function σ(F_1::Rational, m_1::Rational, F_2::Rational, m_2::Rational)
    ψe = Fm_state(F_1, m_1)⊕ Ket(SpinBasis(F_2))
    ψg = Ket(SpinBasis(F_1)) ⊕ Fm_state(F_2, m_2)
    return ψe ⊗ dagger(ψg)
end

"""
Atomic lowering operator for e-g system with hyperfine structure. Spherical basis
"""
function Σ_q(q::Int, F_e::Rational, F_g::Rational)
    # Prevent the case where q-m is larger than F_e
    return dagger(sum(clebschgordan(F_g, -m, 1, q, F_e)*σ(F_e, m - q, F_g, m) for m = -F_g:1:F_g if abs(q-m) <= F_e ))
end