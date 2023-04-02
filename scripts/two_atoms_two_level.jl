using DrWatson
@quickactivate "DecayDynamics"
include("../src/Base.jl")
include("../src/many_atoms_two_level.jl")

# Computing collective effects of two atoms with two-level hyperfine structure

# F_i = [1//1, 0//1]
F_i = [0//1, 1//1]
g_i = [1, 1]
Γ_i = [1, 0]
# x, y, z
Bfield = [0, 0, 0] 
positions = [[0, 0, 0], [0.0, 0, 0.0]] # length of the vector is the number of atoms
# tspan = range(0, -log(0.00000001)/Γ_i[1]/2, length=100)
tspan = range(0, 3, length=100)

ρ_0 = sparse(
    projector((Fm_state(F_i[1], F_i[1])) ⊕ Ket(SpinBasis(F_i[2]))) ⊗ 
    projector((Fm_state(F_i[1], F_i[1])) ⊕ Ket(SpinBasis(F_i[2])))
) # Fully excited state

# ρ_0 = sparse(
#     projector(Fm_state(F_i[1], F_i[1]) ⊕ Ket(SpinBasis(F_i[2]))) ⊗ 
#     projector(Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], F_i[2]))
# ) # ptl1 in excited state


parameters = @dict F_i g_i Γ_i Bfield positions ρ_0 tspan 
@time result = evolve_master(parameters,)

fig = plot_dynamics(result)
plot!(tspan*Γ_i[1], exp.(-Γ_i[1]*tspan), lab=L"\propto e^{-\Gamma t}", color=:black, lw=1)
plot!(tspan*Γ_i[1], exp.(-(2)*Γ_i[1]*tspan), lab=L"\propto e^{-2\Gamma t}", color=:black, ls=:dash, lw=1)
# hline!(fig, [0.1, exp(-tspan[end]*Γ_i[1])], lab="0.1", color=:black, lw=1)
fig