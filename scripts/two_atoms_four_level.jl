using DrWatson
@quickactivate "DecayDynamics"
include("../src/Base.jl")
include("../src/many_atoms_four_level.jl")

# Computing collective effects of two atoms with three level hyperfine structure
F_i = [2//1, 1//1, 1//1, 0//1] # two intermediate states
g_i = [1, 1, 1, 1]
Γ_ij = [1, 0.1, 0] # Reference is 3D1 level.

knorm = 2π*[1, 2.6/0.689, 0] # length of the k-vector, reference to 2.6. index dictate which transition.
Bfield = [0, 0, 0] 
positions = [[100, 0, 0], [0.0, 0, 0]] # 1 = 2.6 um. 
tspan = range(0, 20, length=300)
m_exc = F_i[1] # excited state

# ρ_0 = ( (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3]))) ⊗ 
#         (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3]))) )
# description="m=$(AbstractFloat(m_exc))_FullyExcitedState_r12=$(round.(positions[1]-positions[2], digits=2))"

# ρ_0 = ( (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3]))) ⊗ 
#         ( Ket(SpinBasis(F_i[1])) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Fm_state(F_i[3], F_i[3])) )
# description="m=$(Int(m_exc))_Ptl1OnlyExcied_r12=$(round.(positions[1]-positions[2], digits=2))"

ρ_0 = ( (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3]))) ⊗ 
        (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], F_i[2]) ⊕ Ket(SpinBasis(F_i[3])) ) )
description="m=$(Int(m_exc))_Ptl1ePtl2f_r12=$(round.(positions[1]-positions[2], digits=2))"


parameters = @dict F_i g_i Γ_ij Bfield positions ρ_0 tspan m_exc knorm
@time result = evolve_master(parameters, check_rates=false)
fig = plot_dynamics(result)
filename = savename(tostringdict(parameters), connector="-", allowedtypes=[AbstractFloat, Int, String])
savefig(fig, plotsdir("two_atoms_three_level", description*filename*".pdf"))
