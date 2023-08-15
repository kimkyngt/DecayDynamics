using DrWatson
@quickactivate "DecayDynamics"
include("../src/single_atom_three_level.jl")

# Parameters

F_i = [11//2, 11//2, 9//2]
g_i = [1, 1, 0.0]
Γ_i = [10.0, 1.0, 0]
Bfield = [0, 0, 3]
tspan = 0:0.002:5
frac = normalize([1, 1])

ρ_0 = sparse(normalize(
    projector((sqrt(frac[1])*Fm_state(F_i[1], 5//2)+ sqrt(frac[2])*Fm_state(F_i[1], 3//2)) ⊕ Ket(SpinBasis(F_i[2])) ⊕ Ket(SpinBasis(F_i[3])))  
    ))

parameters = @dict F_i g_i Γ_i Bfield ρ_0 tspan frac

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
    accesses=[:frac, :F_i, :g_i, :Bfield, :Γ_i],
)

result = evolve_master(parameters)
@pack! result = filename

# Save data and gif
# anim = gif( get_animation(result), datadir("three_level_hyperfine_single_atom", filename*".gif"), fps=5)
# savefig(plot_dynamics(result), datadir("three_level_hyperfine_single_atom", filename*".png"))
# wsave(datadir("three_level_hyperfine_single_atom", filename*".jld2"), tostringdict(result))

figs_population = plot_dynamics(result)
figs_detector = plot_detector_signal(result, pi/2-10/180*pi, 0)
figs_detector2 = plot_detector_signal(result, pi/2, 0)
fig = vcat(figs_population, figs_detector, plot(axis=false), figs_detector2)
fig_tot = plot(fig..., size=(1800, 1800), margins=40Plots.px, layout=(4, 3))

# savefig(fig_tot, plotsdir("three_level_hyperfine_single_atom", filename*".pdf"))
display(fig_tot)
