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
positions = [[1/12, 0, 0], [0, 0, 0]] # length of the vector is the number of atoms
# tspan = range(0, -log(0.00000001)/Γ_i[1]/2, length=100)
tspan = range(0, 4, length=500)

m_exc = F_i[1]*0# excited state

ρ_0 = (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2]))) ⊗ (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], F_i[2]))
description="m=$(Int(m_exc))_Atom1onlyExcited_r12=$(round.(positions[1]-positions[2], digits=2))"


parameters = @dict F_i g_i Γ_i Bfield positions ρ_0 tspan m_exc 
@time result = evolve_master(parameters, check_rates=true)
fig = plot_dynamics(result)
filename = savename(tostringdict(parameters), connector="-", allowedtypes=[AbstractFloat, Int, String])
# savefig(fig, plotsdir("Concurrence.pdf"))
fig