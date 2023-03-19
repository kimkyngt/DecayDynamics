using DrWatson
@quickactivate "DecayDynamics"
include("../src/single_atom_three_level.jl")

# Parameters

F_i = [11//2, 9//2, 9//2]
g_i = [127/85, 1, 0.0]
Γ_i = [10.0, 1.0, 0]
Bfield = [0, 0, 2]
tspan = 0:0.02:2

ρ_0 = sparse(normalize(
    projector((Fm_state(F_i[1], 9//2)+ 0*Fm_state(F_i[1], 7//2)+ Fm_state(F_i[1], 5//2)) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3])))  
    ))
file_name_prefix = "test"


parameters = @dict F_i g_i Γ_i Bfield ρ_0 tspan 
filename = file_name_prefix*savename(tostringdict(parameters), connector="-", allowedtypes=[AbstractFloat, Int, String])
result = evolve_master(parameters)
@pack! result = filename

# Save data and gif
# anim = gif( get_animation(result), datadir("three_level_hyperfine_single_atom", filename*".gif"), fps=5)
# savefig(plot_dynamics(result), datadir("three_level_hyperfine_single_atom", filename*".pdf"))
# savefig(plot_dynamics(result), datadir("three_level_hyperfine_single_atom", filename*".png"))
# wsave(datadir("three_level_hyperfine_single_atom", filename*".jld2"), tostringdict(result))

figs_population = plot_dynamics(result)
figs_detector = plot_detector_signal(result, pi/2-0/180*pi, 0)
fig = vcat(figs_population, figs_detector)
plot(fig..., size=(1800, 1500), margins=40Plots.px)