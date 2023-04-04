using DrWatson
@quickactivate "DecayDynamics"
include("../src/Base.jl")
include("../src/many_atoms_four_level.jl")

# Computing collective effects of two atoms with three level hyperfine structure
F_i = [2//1, 1//1, 1//1, 0//1] # two intermediate states, 3D1, 3P1, 3P0, 1S0
g_i = [1, 1, 1, 1]
Γ_kl = sparse(
        [1, 1, 2],
        [2, 3, 4],
        [1*0.385, 1*0.597, 0.1], # Branching ratio
        4, 4
)
knorm = 2π*sparse(
        [1, 1, 2],
        [2, 3, 4],
        [1, 2.6/2.7, 2.6/0.689], # wavevector ratio
        4, 4
) # length of the k-vector, reference to 2.6. index dictate which transition.

Bfield = [0, 0, 0] 
positions = [[0.0, 0, 0], [0.0, 0, 0]] # 1 = 2.6 um. 
tspan = range(0, 10/Γ_kl[1, 2], length=300)
m_exc = F_i[1] # excited state

ρ_0 = ( ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...) ⊗ 
        ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...) )
description="m=$(AbstractFloat(m_exc))_FullyExcitedState_r12=$(round.(positions[1]-positions[2], digits=2))"

# ρ_0 = ( (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3]))) ⊗ 
#         ( Ket(SpinBasis(F_i[1])) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Fm_state(F_i[3], F_i[3])) )
# description="m=$(Int(m_exc))_Ptl1OnlyExcied_r12=$(round.(positions[1]-positions[2], digits=2))"

# ρ_0 = ( (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3]))) ⊗ 
#         (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], F_i[2]) ⊕ Ket(SpinBasis(F_i[3])) ) )
# description="m=$(Int(m_exc))_Ptl1ePtl2f_r12=$(round.(positions[1]-positions[2], digits=2))"


parameters = @dict F_i g_i Γ_kl Bfield positions ρ_0 tspan m_exc knorm
@time result = evolve_master(parameters, check_rates=false)
plot_dynamics(result)
# filename = savename(tostringdict(parameters), connector="-", allowedtypes=[AbstractFloat, Int, String])
# savefig(fig, plotsdir("two_atoms_three_level", description*filename*".pdf"))

# Convet a string to a WignerSymbols