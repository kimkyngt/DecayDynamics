using DrWatson
@quickactivate "DecayDynamics"
include("../src/Base.jl")
include("../src/many_atoms_four_level.jl")

foldername = "two_atoms_four_level"
# subfoldername = "fraction_scan"

# Computing collective effects of two atoms with three level hyperfine structure
F_i = [2//1, 1//1, 1//1, 0//1] # two intermediate states, 3D1, 3P1, 3P0, 1S0
g_i = [1, 1, 1, 1]
Γ_kl = sparse(
        [1, 1, 2],
        [2, 3, 4],
        [1*0.385, 1*0.597, 0.1], # Branching ratio
        4, 4
)
knorm = 2π*sparse(
        [1, 1, 2],
        [2, 3, 4],
        [1, 2.6/2.7, 2.6/0.689], # wavevector ratio
        4, 4
) # length of the k-vector, reference to 2.6. index dictate which transition.
frac = 1
exc_frac = [frac, 1-frac] # Excitation fraction of level 1 and 3
Bfield = [0, 0, 0]
positions = [[0, 0, zz], [0.0, 0, 0]] # 1 = 2.6 um. 
tspan = range(0, 37*0.385/Γ_kl[1, 2], length=300)
m_exc = F_i[1] # excited state

# ρ_0 = ( ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...) ⊗ 
#         ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...) )
# description="[$(file_counter)]m=$(AbstractFloat(m_exc))_FullyExcitedState_r12=$(round.(positions[1]-positions[2], digits=2))"

ψsingle=sqrt(exc_frac[1])*spinup(SpinBasis(F_i[1])) ⊕ Ket(SpinBasis(F_i[2])) ⊕ exp(im*knorm[1, 2]*(norm(positions[1]-positions[2])))*sqrt(exc_frac[2])*spinup(SpinBasis(F_i[3])) ⊕ Ket(SpinBasis(F_i[4]))
ρ_0 = normalize(ψsingle ⊗ ψsingle)
description="[$(file_counter)]m=$(AbstractFloat(m_exc))_VaryingExcitation_r12=$(round.(positions[1]-positions[2], digits=2))"

parameters = @dict F_i g_i Γ_kl Bfield positions ρ_0 tspan m_exc knorm exc_frac
@time result = evolve_master(parameters, check_rates=false)
fig = plot_dynamics(result)
filename = savename(tostringdict(parameters), connector="-", allowedtypes=[AbstractFloat, Int, String, Vector{Int}, Vector{Float64}])
# savefig(fig, plotsdir(foldername, subfoldername, description*filename*".pdf"))
