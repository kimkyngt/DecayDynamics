using DrWatson
@quickactivate "DecayDynamics"
include(srcdir("single_atom_three_level.jl"))

# Parameters

F_i = [11//2, 9//2, 9//2]
g_i = [127/85, 1, 0.0]
Γ_i = [10.0, 1.0, 0]
Bfield = [0, 0, 1]
tspan = 0:0.05:5

# ρ_0 = sparse(normalize(
#     projector(Fm_state(F_i[1], F_i[1]) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3])))  
#     ))
# file_name_prefix = "[Streched]"

# ρ_0 = sparse(normalize(
#     projector(Fm_state(F_i[1], 7//2) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3])))  
#     ))
# file_name_prefix = "[m=7|2]"

# ρ_0 = sparse(normalize(
#     projector((Fm_state(F_i[1], 9//2) + Fm_state(F_i[1], 5//2)) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3])))  
#     ))
# file_name_prefix = "[1m=9|2,1m=5|2]"

# ρ_0 = sparse(normalize(
#     projector((Fm_state(F_i[1], 9//2) + Fm_state(F_i[1], 7//2)) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3])))  
#     ))
# file_name_prefix = "[1m=9|2,1m=7|2]"

# ρ_0 = sparse(normalize(
#     projector((Fm_state(F_i[1], 9//2) + Fm_state(F_i[1], -9//2)) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3])))  
#     ))
# file_name_prefix = "[1m=9|2,1m=-9|2]"

ρ_0 = sparse(normalize(
    projector((Fm_state(F_i[1], 9//2) + Fm_state(F_i[1], 7//2) + Fm_state(F_i[1], 5//2)) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3])))  
    ))
file_name_prefix = "[1m=9|2,1m=7|2,1m=5|2]"


parameters = @dict F_i g_i Γ_i Bfield ρ_0 tspan 
filename = file_name_prefix*savename(tostringdict(parameters), connector="-", allowedtypes=[AbstractFloat, Int, String])
@time result = evolve_master(parameters)
@pack! result = filename

# Save data and gif
# anim = gif( get_animation(result), datadir("three_level_hyperfine_single_atom", "examples", filename*".gif"), fps=5)
# savefig(plot_dynamics(result), datadir("three_level_hyperfine_single_atom", "examples", filename*".pdf"))
# savefig(plot_dynamics(result), datadir("three_level_hyperfine_single_atom", "examples", filename*".png"))
# wsave(datadir("three_level_hyperfine_single_atom", "examples", filename*".jld2"), tostringdict(result))

plot_dynamics(result)