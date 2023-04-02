using DrWatson
@quickactivate "DecayDynamics"
include("../src/Base.jl")
include("../src/many_atoms_two_level.jl")

# Computing collective effects of two atoms with two-level hyperfine structure

# F_i = [1//1, 0//1]
F_i = [1//1, 0//1]
g_i = [1, 1]
Γ_i = [1, 0]
# x, y, z
Bfield = [0, 0, 0] 
positions = [[0, 0, 0], [0, 0, 0.15]] # length of the vector is the number of atoms
# tspan = range(0, -log(0.00000001)/Γ_i[1]/2, length=100)
tspan = range(0, 10, length=500)

m_exc = F_i[1]-1 # excited state

# ρ_0 = (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2]))) ⊗ (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2]))) 
# description="Fully-excited-state_r12=$(round.(positions[1]-positions[2], digits=2))"

ρ_0 = (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2]))) ⊗ (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], F_i[2]))
description="Atom1onlyExcited_r12=$(round.(positions[1]-positions[2], digits=2))"

parameters = @dict F_i g_i Γ_i Bfield positions ρ_0 tspan m_exc 
@time result = evolve_master(parameters, check_rates=true)
fig = plot_dynamics(result)
filename = savename(tostringdict(parameters), connector="-", allowedtypes=[AbstractFloat, Int, String])
savefig(fig, plotsdir("two_atoms_two_level", description*filename*".pdf"))

fig