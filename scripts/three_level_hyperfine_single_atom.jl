using DrWatson
@quickactivate "DecayDynamics"
include("../src/single_atom_three_level.jl")

# Parameters
filename = "test"

F_i = [3//2, 3//2, 3//2]
g_i = [1.0, 1.0, 1.0]
Γ_i = [10.0, 1.0, 0]
Bfield = [0, 0, 1]
tspan = 0:0.05:3
ρ_0 = sparse(projector(normalize(Fm_state(F_i[1], F_i[1]) + Fm_state(F_i[1], F_i[1]-1) - Fm_state(F_i[1], F_i[1]-2)) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3]))))

parameters = @dict F_i g_i Γ_i Bfield ρ_0 tspan 

@time result = evolve_master(parameters)
@pack! result = filename

plot_dynamics(result)

# Save data and gif
anim = gif( get_animation(result), datadir("three_level_hyperfine_single_atom", filename*".gif"), fps=5)
savefig(plot_dynamics(result), datadir("three_level_hyperfine_single_atom", filename*".pdf"))
wsave(datadir("three_level_hyperfine_single_atom", filename*".jld2"), tostringdict(result))

