using DrWatson, DataFrames
import Rotations: RotX, RotZ
@quickactivate "DecayDynamics"
include("../../src/Base.jl")
include("../../src/many_atoms_four_level.jl")

subfoldername = "m=0.5_fraction_scan_0.951S0_r12=2.6um"
DATA_DIR = datadir("two_atoms_four_level", subfoldername)
PLOTS_DIR = plotsdir("two_atoms_four_level", subfoldername)

df = collect_results(DATA_DIR)

function get_pop_int_comparison(r_det, fileindx)
    println("Processing file $(df[fileindx, :path])")

    # Get parameters
    Γ3D1 = df[fileindx, "Γ_kl"][1, 3] + df[fileindx, "Γ_kl"][1, 2]
    ρ_t = df[fileindx, "ρ_t"]
    F_i = df[fileindx, "F_i"]
    tΓ = df[fileindx, "t_out"]*Γ3D1
    k689 = df[1, "knorm"][2, 4]
    positions = df[fileindx, "positions"]

    # compute 3P1 population
    pop3P1 = [real.(expect(  
                ⊕([kk == 2 ? identityoperator(SpinBasis(F_i[kk])) : projector(Ket(SpinBasis(F_i[kk]))) for kk in eachindex(F_i)]...) ⊗ ⊕([identityoperator(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
                , ρ_t[jj]) ) for jj in eachindex(ρ_t)
    ]
    ypop = normalize(pop3P1)
    fig1 = plot(xlab="tΓ", ylab="Area-normalized signal")
    fig2 = plot(xlab="tΓ", ylab="population - intensity", legend=false)

    # compute intensity

    θarr = range(0, π/2, length=11)
    ϕarr = 0
    mycolor = cgrad(:roma, length(θarr)*length(ϕarr), categorical=true)
    colorindx = 1
    for θ in θarr
        for ϕ in ϕarr
            r_det_pos = r_det*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
                # compute intensity                     
            Intensity = [real.(sum(ϵ_q(p)*(GreenTensor(r_det_pos - positions[i], k689)'*GreenTensor(r_det_pos - positions[j], k689)) * adjoint(ϵ_q(q)) * expect(dagger(Σ_iq(i, p, F_i, 2, 4))*Σ_iq(j, q, F_i, 2, 4), ρ_t[timeindx]) for p in [-1, 0, 1] for q in [-1, 0, 1] for i in [1, 2] for j in [1, 2])) for timeindx in eachindex(ρ_t)]

            yint = normalize(Intensity)
            plot!(fig1, tΓ, yint, lab="θ=$(round(θ/π, digits=2))π, ϕ=$(round(ϕ/π, digits=2))π", color=mycolor[colorindx])
            plot!(fig2, tΓ, ypop - yint, lab="θ=$(round(θ/π, digits=2))π, ϕ=$(round(ϕ/π, digits=2))π", color=mycolor[colorindx])
            colorindx += 1
        end
    end

    plot!(fig1, tΓ, ypop, lab="population", ls=:dash, color=:black)
    plot(fig1, fig2, layout=(2, 1), size=(800, 800), plot_title="3D1-3P0-1S0 pop=$(round.([df[fileindx, :exc_frac]; 1-sum(df[fileindx, :exc_frac])], digits=3))")
end


# Analysis parameters
r_det = 10 # position of detector

for fileindx in axes(df, 1)
    fig = get_pop_int_comparison(r_det, fileindx)
    savefig(fig, joinpath(PLOTS_DIR, splitdir(df[fileindx, :path])[2]*"_int-pop.pdf"))
    display(fig)
end


# savefig(fig, joinpath(PLOTS_DIR, splitdir(df[fileindx, :path])[2]*"_radiation_pattern_logInt.pdf"))