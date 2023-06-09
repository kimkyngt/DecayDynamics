using DrWatson
@quickactivate "DecayDynamics"
include(srcdir("single_atom_two_level.jl"))

# Parameters
F_e = 1//1
F_g = 0//1
tspan = 0:0.05:5
# Bfield in units of Γ/2π
field_x = 0.0
field_y = 0.0
field_z = 1
# g-factor
g_e = 1.0 
g_g = 1.0

# state
# ρ_0 = sparse(projector(normalize( Fm_state(F_e, F_e)) ⊕ Ket(SpinBasis(F_g))))
# file_name_prefix="[Example]m_e=1"

# ρ_0 = sparse(projector(normalize( Fm_state(F_e, F_e-1)) ⊕ Ket(SpinBasis(F_g))))
# file_name_prefix="[Example]m_e=0"

# ρ_0 = sparse(projector(normalize( Fm_state(F_e, F_e) + Fm_state(F_e, F_e-1)) ⊕ Ket(SpinBasis(F_g))))
# file_name_prefix="[Example]m_e=0|m_e=1|"

# ρ_0 = sparse(normalize( projector(Fm_state(F_e, F_e) ⊕ Ket(SpinBasis(F_g))) + projector(Fm_state(F_e, F_e-1) ⊕ Ket(SpinBasis(F_g))) ))
# file_name_prefix="[Example]mixed|m_e=0|m_e=1|"

# ρ_0 = sparse(normalize( projector(Fm_state(F_e, F_e) ⊕ Ket(SpinBasis(F_g))) + 0.5*projector(Fm_state(F_e, F_e-1) ⊕ Ket(SpinBasis(F_g))) ))
# file_name_prefix="[Example]mixed|1m_e=0|0.5m_e=1|"

# ρ_0 = sparse(projector(normalize( Fm_state(F_e, F_e) + Fm_state(F_e, F_e-2)) ⊕ Ket(SpinBasis(F_g))))
# file_name_prefix="[Example]m_e=-1|m_e=1|"

# ρ_0 = sparse(projector(normalize( Fm_state(F_e, F_e)  + Fm_state(F_e, F_e-1) - Fm_state(F_e, F_e-2)) ⊕ Ket(SpinBasis(F_g))))
# file_name_prefix="[Example]m_e=-1m_e=0|-1m_e=1|"

ρ_0 = sparse(projector(normalize( Fm_state(F_e, F_e)  + 2*Fm_state(F_e, F_e-1) - Fm_state(F_e, F_e-2)) ⊕ Ket(SpinBasis(F_g))))
file_name_prefix="[Example]m_e=-1|2m_e=0|-1m_e=1|"

parameters = @dict F_e F_g ρ_0 tspan field_x field_z field_y g_e g_g 
@time result = evolve_master(parameters)
filename = file_name_prefix*savename(tostringdict(parameters), connector="-", allowedtypes=[AbstractFloat, Int, String])
@pack! result = filename

fig = plot_dynamics(result)

# Save data and gif
anim = gif( get_animation(result), datadir("two_level_hyperfine_single_atom", "examples", filename*".gif"), fps=5)
savefig(fig, datadir("two_level_hyperfine_single_atom", "examples", filename*".pdf"))
savefig(fig, datadir("two_level_hyperfine_single_atom", "examples", filename*".png"))
wsave(datadir("two_level_hyperfine_single_atom", "examples", filename*".jld2"), tostringdict(result))

