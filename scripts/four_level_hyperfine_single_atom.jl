# initialize the project
using DrWatson
@quickactivate "DecayDynamics"
include("../src/Base.jl")
include("../src/single_atom_four_level.jl")
# Simulation parameters
F_i = [1//2, 1//2, 1//2, 1//2] # two intermediate states, 3D1, 3P1, 3P0, 1S0
g_i = [1, 1, 1, 1]
Γ_kl = sparse(
        [1, 1, 2],
        [2, 3, 4],
        [0*0.385, 0*0.597, 0.1], # Branching ratio, 3P1 assumed to be just 10 times slower. 
        4, 4
)
frac = 1
exc_frac = [frac, 1-frac] # Excitation fraction of level 1 and 3
Bfield = [0, 0, 0.1]
tspan = range(0, 20, length=300)
ρ_0 = normalize(directsum(
    sqrt(exc_frac[1])*(Fm_state(F_i[1], 1//2) + Fm_state(F_i[1], -1//2)) ,
    Ket(SpinBasis(F_i[2])),
    sqrt(exc_frac[2])*spinup(SpinBasis(F_i[3])),
    Ket(SpinBasis(F_i[4]))
))

parameters = @dict F_i g_i Bfield Γ_kl ρ_0 tspan exc_frac

result = evolve_master(parameters)

figs_population = plot_dynamics(result)
plot(figs_population..., size=(2000, 1200), layout=(3, 3), margins=10Plots.mm)

