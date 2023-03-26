using DrWatson
@quickactivate "DecayDynamics"
include("../src/single_atom_two_level.jl")

# Parameters
F_e = 1//1
F_g = 0//1
Fe = string(F_e.num)*"half"
Fg = string(F_g.num)*"half"
ρ_0 = sparse(projector(normalize( Fm_state(F_e, F_e) + 1*Fm_state(F_e, F_e-1)) ⊕ Ket(SpinBasis(F_g))))
# ρ_0 = sparse(normalize(one(SpinBasis(F_e))) ⊕ projector(Ket(SpinBasis(F_g))))
# ρ_0 = sparse(normalize(projector(Fm_state(F_e, F_e)⊕ Ket(SpinBasis(F_g))) + projector(Fm_state(F_e, F_e-1)⊕ Ket(SpinBasis(F_g)))))
Astring="1_3_superposition2_gratio3"
tspan = 0:0.2:5
# Bfield in units of Γ/2π
field_x = 0.0
field_y = 0.0
field_z = 1
# g-factor
g_e = 1.0 
g_g = 0.3

parameters = @dict F_e F_g Fe Fg ρ_0 tspan field_x field_z field_y g_e g_g Astring 
@time result = evolve_master(parameters)
# @time draw_radiation_pattern(ρ_0, title="Initial state", size=(400, 400), camera=[45, 25])
filename = savename(tostringdict(parameters), connector="-", allowedtypes=[AbstractFloat, Int, String])
@pack! result = filename


# Save data and gif
# gr()
# anim = gif( get_animation(result), datadir("two_level_hyperfine_single_atom", filename*".gif"), fps=5)
# savefig(plot_dynamics(result), datadir("two_level_hyperfine_single_atom", filename*".pdf"))
# wsave(datadir("two_level_hyperfine_single_atom", filename*".jld2"), tostringdict(result))
gr()
plot_dynamics(result)
gif( get_animation(result), fps=5)