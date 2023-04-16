using DrWatson, DataFrames
import Rotations: RotX, RotZ
@quickactivate "DecayDynamics"
include("../../src/Base.jl")
include("../../src/many_atoms_four_level.jl")

subfoldername = "F_i=[2,1,1,0]position_scan_wide"
DATA_DIR = datadir("two_atoms_four_level", subfoldername)

# Calculate Rabi frequency ratios 
if ~(@isdefined df)
    df = collect_results(DATA_DIR)
else
    print("df exist, skip loading \n")
end

function main(xyz::String, m_exc::Int64)
    println("Plotting m=$(m_exc), $(xyz)-direction scan")
    df_m2 = df[df[:, "m_exc"] .== m_exc, :]
    df_m2[:, "direction"] =[df_m2[ii, "displacement"][3] == 0 ? "x" : "z" for ii in range(1, nrow(df_m2))]
    df_m2[:, "rij"] =[norm(df_m2[ii, "displacement"]) for ii in range(1, nrow(df_m2))]
    df_m2_z = df_m2[df_m2[:, "direction"] .== xyz, :]
    sort!(df_m2_z, ["rij"])

    # plot observables
    timeindice = [2, 5, 10, 20, 60,100]
    figs = []
    fileindx = 10
    @showprogress for timeindx in timeindice
        Γ3D1 = df_m2_z[timeindx, "Γ_kl"][1, 3] + df_m2_z[timeindx, "Γ_kl"][1, 2]
        # Get parameters
        k689 = df_m2_z[1, "knorm"][2, 4]
        positions = df_m2_z[fileindx, "positions"]
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
                Intensity_xy[xindx, yindx] = real.(sum(ϵ_q(p)*(GreenTensor(r_det_xy - positions[i], k689)'*GreenTensor(r_det_xy - positions[j], k689)) * adjoint(ϵ_q(q)) * expect(dagger(Σ_iq(i, p, df_m2_z[fileindx, "F_i"], 2, 4))*Σ_iq(j, q, df_m2_z[fileindx, "F_i"], 2, 4), df_m2_z[fileindx, "ρ_t"][timeindx]) for p in [-1, 0, 1] for q in [-1, 0, 1] for i in [1, 2] for j in [1, 2]))
            end
        end

        # compute in xz-plane
        Intensity_xz = zeros(Nx, Nz)
        for xindx in eachindex(x_list)
            for zindx in eachindex(z_list)
                r_det_xz = [x_list[xindx], 0, z_list[zindx]]
                Intensity_xz[xindx, zindx] = real.(sum(ϵ_q(p)*(GreenTensor(r_det_xz - positions[i], k689)'*GreenTensor(r_det_xz - positions[j], k689)) * adjoint(ϵ_q(q)) * expect(dagger(Σ_iq(i, p, df_m2_z[fileindx, "F_i"], 2, 4))*Σ_iq(j, q, df_m2_z[fileindx, "F_i"], 2, 4), df_m2_z[fileindx, "ρ_t"][timeindx]) for p in [-1, 0, 1] for q in [-1, 0, 1] for i in [1, 2] for j in [1, 2]))
            end
        end
        fig1 =  contour(x_list, y_list, log10.(Intensity_xy'), title="tΓᴰ=$(round(df_m2_z[fileindx, "t_out"][timeindx]*Γ3D1, sigdigits=2)) r₁₂=$(round(norm(positions[1]-positions[2]), sigdigits=2)), xy-plane", margins=10Plots.mm, aspect_ratio=1, legend=false, xlab="x/(2.6 μm)", ylab="y/(2.6 μm)", fill=true, levels=10, lw=0, color=reverse(cgrad(:oslo)))
        scatter!(fig1, [positions[1][1]], [positions[1][2]], label="atom 1", markersize=5, color=:white)
        scatter!(fig1, [positions[2][1]], [positions[2][2]], label="atom 2", markersize=5, color=:white)

        fig2 =  contour(x_list, z_list, log10.(Intensity_xz'), title="xz-plane", margins=10Plots.mm, aspect_ratio=1, legend=false, xlab="x/(2.6 μm)", ylab="z/(2.6 μm)", fill=true, levels=10, lw=0, color=reverse(cgrad(:oslo)))
        scatter!(fig2, [positions[1][1]], [positions[1][3]], label="atom 1", markersize=5, color=:white)
        scatter!(fig2, [positions[2][1]], [positions[2][3]], label="atom 2", markersize=5, color=:white)
        
        push!(figs, plot(fig1, fig2, layout=(1, 2), size=(800, 400), ))
    end
    return plot(figs..., layout=(3, 2), size=(1600, 1200))
end

main("z", 1)