using DrWatson
@quickactivate "DecayDynamics"
include("../src/base.jl")
include("../src/many_atoms_two_level.jl")

# Computing collective effects of two atoms with two-level hyperfine structure

F_i = [3//2, 1//2]
g_i = [1, 1]
Γ_i = [1, 0]
# x, y, z
Bfield = [0, 0, 1] 
kr = [[0, 0, 0], [0, 0, 1]] # length of the vector is the number of atoms
tspan = 0:0.02:2

ρ_0 = sparse(normalize(
    projector((Fm_state(F_i[1], F_i[1])) ⊕ Ket(SpinBasis(F_i[2]))) ⊗ 
    projector((Fm_state(F_i[1], F_i[1])) ⊕ Ket(SpinBasis(F_i[2])))
)) # Fully excited state

parameters = @dict F_i g_i Γ_i Bfield kr ρ_0 tspan 
@time result = evolve_master(parameters)

fig = plot_dynamics(result)
fig