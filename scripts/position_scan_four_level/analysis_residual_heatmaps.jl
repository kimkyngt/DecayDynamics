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


function get_idir_plot(xyz::String, m_exc::Int64)
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

        # Plot raw data
        # plot!(fig1, xdata, ydata, alpha=0.1, color=:black, lw=1)

        # Fit data
        @. model(x, p) = p[1]*(exp(-(x - p[4])/p[2]) - exp(-(x - p[4])/p[3])) + p[5]
        transition_indx = findall(!iszero, Γ_kl)
        p0 = [maximum(ydata), 1/Γ_kl[transition_indx[3]], τtot, 0, 0] # True taus
        fit_raw = curve_fit(model, xdata, ydata, p0)
        fit_model = model(xdata, fit_raw.param)

        residual = ydata - fit_model
        # plot!(fig2, tΓ, residual, label="data - fit", alpha=0.1, color=:black, lw=1)
        residuals[ii, :] = residual

    end
    fig = plot(
        time, distance, residuals*1e3, 
        st=:heatmap,  
        clims=(-1, 1), 
        color=:vik,
        fill=true,
        xlabel="tΓ", 
        ylab="$(xyz)-distance/(2.6 μm)", 
        title="m=$(Int64(m_exc)), Fit residuals × 10³", 
        yticks=0:0.25:1,
        xticks=[0, 2, 4, 6, 8, 10, 20, 30], 
        grid=true
        )
    return fig
end

fig1 = plot(get_idir_plot("x", 2), get_idir_plot("z", 2), layout=(2, 1), size=(800, 800))
display(fig1)
savefig(fig1, plotsdir("m=2,"*subfoldername*".pdf"))
savefig(fig1, plotsdir("m=2,"*subfoldername*".png"))


fig2 = plot(get_idir_plot("x", 1), get_idir_plot("z", 1), layout=(2, 1), size=(800, 800))
savefig(fig2, plotsdir("m=1,"*subfoldername*".pdf"))
savefig(fig2, plotsdir("m=1,"*subfoldername*".png"))

