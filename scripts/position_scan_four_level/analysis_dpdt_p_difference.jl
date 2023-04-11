using DrWatson, DataFrames
@quickactivate "DecayDynamics"
include("../../src/Base.jl")
include("../../src/many_atoms_four_level.jl")

subfoldername = "F_i=[2,1,1,0]position_scan"
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
    sort!(df_m2_z, [:rij])

    # plot observables
    time = undef
    distance = []
    residuals = Matrix(undef, size(df_m2_z)[1], length(df_m2_z[1, "t_out"]))
    mycolor = cgrad(:viridis, Int(size(df_m2_z)[1]), categorical=true)
    fig = plot(legend=false)
    for ii in range(1, size(df_m2_z)[1])
        # global fig1, fig2, distance, residuals, time
        # Get parameters
        F_i = df_m2_z[ii, "F_i"]
        ρ_t = df_m2_z[ii, "ρ_t"]
        t_out = df_m2_z[ii, "t_out"]
        Γ_kl = df_m2_z[ii, "Γ_kl"]
        push!(distance, df_m2_z[ii, "rij"])

        # Compute taus
        τtot = 1/(Γ_kl[1, 3] + Γ_kl[1, 2])
        τ1 = 1/Γ_kl[2, 4]
        τ2 = 1/Γ_kl[1, 2]
        tΓ = t_out/τtot#*Γ_kl[1, 2]
        time = tΓ

        # Get 3P1 population after partial trace of the other atom.
        xdata = t_out
        ydata = [real.(expect(⊕([kk == 2 ? identityoperator(SpinBasis(F_i[kk])) : projector(Ket(SpinBasis(F_i[kk]))) for kk in eachindex(F_i)]...), ptrace(ρ_t[jj], 2)) ) for jj in eachindex(ρ_t)]

        # Compute derivative 
        dydt = [ydata[ii] - ydata[ii-1] for ii in range(2, length(ydata))]/(xdata[2] - xdata[1])
        xdata_dydt =         [xdata[ii] for ii in range(2, length(xdata))]
        plot!(fig, xdata, ydata, color=mycolor[ii])
        plot!(fig, xdata_dydt, dydt, color=mycolor[ii], ls=:dash)
    end
    
    return fig
end

main("x", 2)