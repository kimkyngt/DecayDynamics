using DrWatson, Plots
@quickactivate "DecayDynamics"
include("../src/base.jl")

F_i = [5//2, 7//2,]
m1 = 3//2# F_i[1] - 1
ψ = Fm_state(F_i[1], m1) ⊕ Ket(SpinBasis(F_i[2]))
ψ_1 = (Σ_q(-1, F_i, 1, 2) + Σ_q(0, F_i, 1, 2) + Σ_q(1, F_i, 1, 2) ) * ψ 


# ΣΣ = (2*F_i[1] + 1)*sum(dagger(Σ_q(q, F_i, 1, 2))*Σ_q(q, F_i, 1, 2) for q in [-1, 0, 1])

ΣΣ = sum(dagger(Σ_q(q, F_i, 1, 2))*Σ_q(q, F_i, 1, 2) for q in [-1, 0, 1])
 