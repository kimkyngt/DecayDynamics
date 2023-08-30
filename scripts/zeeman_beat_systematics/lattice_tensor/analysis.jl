using DrWatson, DataFrames, Statistics
include(srcdir("Base.jl"))
include(srcdir("single_atom_four_level.jl"))

# load data if `df` not exists
if !isdefined(Main, :df)
    df = collect_results(datadir("zeeman_beat", "lattice_tensor"))
end


function get_fit(df, indx)
    df_result = df[indx, :result]
    It_avg = df_result["It"]
    t_out = df_result["t_out"]
    Bfield = df[indx, :params]["Bfield"]
    thetaDet = df[indx, :params]["thetaDet"]
    # fit the averaged data
    @. model(x, p) = p[1]*(exp(-(x - p[4])/p[2]) - exp(-(x - p[4])/p[3])) + p[5] 
    p0 = [maximum(It_avg), 10, 1, 0, 0]
    fit_raw = curve_fit(model, t_out, It_avg, p0)
    fit_model = model(t_out, fit_raw.param)
    se = LsqFit.stderror(fit_raw)

    dtauP = (fit_raw.param[2] - 10)/10
    dtauD = (fit_raw.param[3] - 1)

    fig = plot(t_out, It_avg, lab="Inteisty at θ=$(round(thetaDet/π*180, digits=2))°", dpi=200)
    plot!(t_out, fit_model, lab="double exponential fit", ls=:dash)
    fig_residual = plot(t_out, fit_raw.resid, lab="fit residuals")
    figtot = plot(fig, fig_residual, layout=(2, 1), size=(800, 600), suptitle="ΔτD/τD = $(round(dtauD, digits=4)), ΔτP/τP = $(round(dtauP, digits=3)), B=$(round.(Bfield, digits=2))", titlefontsize=8, xlab="tΓ")
    savefig(figtot, datadir("zeeman_beat", "lattice_tensor", "plot_$(savename(df[indx, :params])).pdf"))
    savefig(figtot, datadir("zeeman_beat", "lattice_tensor", "plot_$(savename(df[indx, :params])).png"))

    return figtot
end

get_fit(df, 1)
get_fit(df, 2)