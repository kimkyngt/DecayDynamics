using DrWatson, DataFrames
import Rotations: RotX, RotZ
@quickactivate "DecayDynamics"
include("../../src/Base.jl")
include("../../src/many_atoms_four_level.jl")

subfoldername = "x-displacement_m=2.5_fraction_scan_0.951S0_r12=10.0um"
DATA_DIR = datadir("two_atoms_four_level", subfoldername)
PLOTS_DIR = plotsdir("two_atoms_four_level", subfoldername)
df = collect_results(DATA_DIR)    

for fileindx in axes(df)[1][5]
    println("Processing file $(df[fileindx, :path])")
    figs = []
    # plot radiation pattern at certain time. 
    timeindice = [2, 10, 60, 300]
    @showprogress for timeindx in timeindice
        Γ3D1 = df[fileindx, "Γ_kl"][1, 3] + df[fileindx, "Γ_kl"][1, 2]
        # Get parameters
        k689 = df[1, "knorm"][2, 4]
        positions = df[fileindx, "positions"]
        Nx = 30
        Ny = 30
        Nz = 30
        x_list = range(-10, 10, length=Nx)
        y_list = range(-10, 10, length=Ny)
        z_list = range(-10, 10, length=Nz)
        Intensity_xy = zeros(Nx, Ny)

        # compute in xy-plane
        for xindx in eachindex(x_list)
            for yindx in eachindex(y_list)
                r_det_xy = [x_list[xindx], y_list[yindx], 0]
                Intensity_xy[xindx, yindx] = real.(sum(ϵ_q(p)*(GreenTensor(r_det_xy - positions[i], k689)'*GreenTensor(r_det_xy - positions[j], k689)) * adjoint(ϵ_q(q)) * expect(dagger(Σ_iq(i, p, df[fileindx, "F_i"], 2, 4))*Σ_iq(j, q, df[fileindx, "F_i"], 2, 4), df[fileindx, "ρ_t"][timeindx]) for p in [-1, 0, 1] for q in [-1, 0, 1] for i in [1, 2] for j in [1, 2]))
            end
        end

        # compute in xz-plane
        Intensity_xz = zeros(Nx, Nz)
        for xindx in eachindex(x_list)
            for zindx in eachindex(z_list)
                r_det_xz = [x_list[xindx], 0, z_list[zindx]]
                Intensity_xz[xindx, zindx] = real.(sum(ϵ_q(p)*(GreenTensor(r_det_xz - positions[i], k689)'*GreenTensor(r_det_xz - positions[j], k689)) * adjoint(ϵ_q(q)) * expect(dagger(Σ_iq(i, p, df[fileindx, "F_i"], 2, 4))*Σ_iq(j, q, df[fileindx, "F_i"], 2, 4), df[fileindx, "ρ_t"][timeindx]) for p in [-1, 0, 1] for q in [-1, 0, 1] for i in [1, 2] for j in [1, 2]))
            end
        end
        fig1 =  contourf(x_list, y_list, log10.(Intensity_xy'), title="tΓᴰ=$(round(df[fileindx, "t_out"][timeindx]*Γ3D1, sigdigits=2)) r₁₂=$(round(norm(positions[1]-positions[2]), sigdigits=2)), xy-plane", aspect_ratio=1, legend=false, xlab="x/(2.6 μm)", ylab="y/(2.6 μm)", levels=10, left_margin=30Plots.px)
        scatter!(fig1, [positions[1][1]], [positions[1][2]], label="atom 1", markersize=5, color=:white)
        scatter!(fig1, [positions[2][1]], [positions[2][2]], label="atom 2", markersize=5, color=:white)

        fig2 =  contourf(x_list, z_list, log10.(Intensity_xz'), title="xz-plane", aspect_ratio=1, legend=false, xlab="x/(2.6 μm)", ylab="z/(2.6 μm)", levels=10, left_margin=30Plots.px)
        scatter!(fig2, [positions[1][1]], [positions[1][3]], label="atom 1", markersize=5, color=:white)
        scatter!(fig2, [positions[2][1]], [positions[2][3]], label="atom 2", markersize=5, color=:white)
        
        push!(figs, plot(fig1, fig2, layout=(1, 2), size=(1000, 400), ))
    end
    fig = plot(figs..., layout=(length(timeindice), 1), size=(1000, 400*length(timeindice)))
    display(fig)
    savefig(fig, joinpath(PLOTS_DIR, splitdir(df[fileindx, :path])[2]*"_radiation_pattern_logInt.pdf"))
end