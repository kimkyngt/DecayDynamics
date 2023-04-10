using DrWatson
@quickactivate "DecayDynamics"
include("../src/Base.jl")
include("../src/many_atoms_four_level.jl")

foldername = "two_atoms_four_level"
subfoldername = "position_scan"

# Computing collective effects of two atoms with three level hyperfine structure
F_i = [5//2, 3//2, 3//2, 1//2] # two intermediate states, 3D1, 3P1, 3P0, 1S0

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
# Bfield = [0, 0, 0]
# g_i = [1, 1, 1, 1]
positions = [[00.02, 0.05, 0.03], [0.0, 0, 0]] # 1 = 2.6 um. 
displacement = positions[1]-positions[2] # for saving
tspan = range(0, 20, length=300)
m_exc = F_i[1]-1  # excited state

# ρ_0 = ( ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...) ⊗ 
#         ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...) )
# description="m=$(AbstractFloat(m_exc))_FullyExcitedState_r12=$(round.(positions[1]-positions[2], digits=2))"

# ψsingle_exc=sqrt(exc_frac[1])*spinup(SpinBasis(F_i[1])) ⊕ Ket(SpinBasis(F_i[2])) ⊕ exp(im*knorm[1, 2]*(norm(positions[1]-positions[2])))*sqrt(exc_frac[2])*spinup(SpinBasis(F_i[3])) ⊕ Ket(SpinBasis(F_i[4]))

ψsingle_exc = ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
ρ_0 = normalize(ψsingle_exc ⊗ ψsingle_exc)
description="Fully_excited"

parameters = @dict F_i Γ_kl positions displacement ρ_0 tspan m_exc knorm exc_frac description
@time result = evolve_master(parameters, check_rates=false)
fig = plot_dynamics(result)
function strfloat(x)
        a = string.(round.(float.(x), digits=3))
        if a isa String
                return a
        else
                return join(a, ",")
        end
end
filename = savename(
        parameters,
        connector="|", 
        allowedtypes=[Real, Int, Vector{<:Number}, ],
        sort=true,
        val_to_string=strfloat,
        accesses=[:description, :F_i, :m_exc, :displacement, :exc_frac],
)

# savefig(fig, plotsdir(foldername, subfoldername, description*filename*".pdf"))
println(filename)
fig