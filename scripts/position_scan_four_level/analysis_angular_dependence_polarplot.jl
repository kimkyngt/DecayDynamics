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
    # fileindice = range(1, size(df_m2_z)[1], step=20)#1:7:60
    fileindice = 1:2:10
    mycolor = cgrad(:roma, fileindice[end], categorical=true)
    fig = plot(legend=true)
    fig1 = plot(legend=true)
    fig2 = plot(legend=true)
    fig3 = plot(legend=true)
    @showprogress for ii in fileindice#
        # Get parameters
        ρ_t = df_m2_z[ii, "ρ_t"]
        t_out = df_m2_z[ii, "t_out"]
        Γ_kl = df_m2_z[ii, "Γ_kl"]
        k689 = df_m2_z[1, "knorm"][2, 4]
        positions = df_m2_z[ii, "positions"]
        # r_det = [1, 2, 2] # detector position
        ϕ_list = range(0, 2π, length=200)
        Intensity_xy = []
        Intensity_xz = []
        Intensity_yz = []
        # for ρ in df_m2_z[ii, "ρ_t"]
        for ϕ in ϕ_list
            r_det = 5
            x, y, z = polar_2_xyz(π/2, ϕ, r_det)
            r_det_xy = [x, y, z] 
            r_det_xz = Vector(RotX(π/2)*[x, y, z]) 
            r_det_yz = Vector(RotZ(π/2)*RotX(π/2)*[x, y, z]) 
            timeindx = 10

            # compute in xy-plane
            push!(Intensity_xy,
                sum(ϵ_q(p)*(GreenTensor(r_det_xy - positions[i], k689)'*GreenTensor(r_det_xy - positions[j], k689)) * adjoint(ϵ_q(q)) * expect(dagger(Σ_iq(i, p, df_m2_z[ii, "F_i"], 2, 4))*Σ_iq(j, q, df_m2_z[ii, "F_i"], 2, 4), df_m2_z[ii, "ρ_t"][timeindx]) for p in [-1, 0, 1] for q in [-1, 0, 1] for i in [1, 2] for j in [1, 2])
            )
            
            # compute in xz-plane
            push!(Intensity_xz,
                sum(ϵ_q(p)*(GreenTensor(r_det_xz - positions[i], k689)'*GreenTensor(r_det_xz - positions[j], k689)) * adjoint(ϵ_q(q)) * expect(dagger(Σ_iq(i, p, df_m2_z[ii, "F_i"], 2, 4))*Σ_iq(j, q, df_m2_z[ii, "F_i"], 2, 4), df_m2_z[ii, "ρ_t"][timeindx]) for p in [-1, 0, 1] for q in [-1, 0, 1] for i in [1, 2] for j in [1, 2])
            )

            push!(Intensity_yz,
                sum(ϵ_q(p)*(GreenTensor(r_det_yz - positions[i], k689)'*GreenTensor(r_det_yz - positions[j], k689)) * adjoint(ϵ_q(q)) * expect(dagger(Σ_iq(i, p, df_m2_z[ii, "F_i"], 2, 4))*Σ_iq(j, q, df_m2_z[ii, "F_i"], 2, 4), df_m2_z[ii, "ρ_t"][timeindx]) for p in [-1, 0, 1] for q in [-1, 0, 1] for i in [1, 2] for j in [1, 2])
            )

        end

        plot!(fig1, ϕ_list, real.(Intensity_xy), proj=:polar, color=mycolor[ii], label=df_m2_z[ii, "rij"], title="xy-plane", legend=:bottomright)
        plot!(fig2, ϕ_list, real.(Intensity_xz), proj=:polar, color=mycolor[ii], label=df_m2_z[ii, "rij"], title="xz-plane", legend=:bottomright)
        plot!(fig3, ϕ_list, real.(Intensity_yz), proj=:polar, color=mycolor[ii], label=df_m2_z[ii, "rij"], title="yz-plane", legend=:bottomright)
    end
    fig = plot(fig1, fig2, fig3, layout=(1, 3), size=(1200, 400))
    return fig
end

main("x", 2)