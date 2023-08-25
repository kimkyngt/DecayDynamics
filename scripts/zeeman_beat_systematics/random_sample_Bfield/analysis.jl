using DrWatson, DataFrames, Statistics
include(srcdir("Base.jl"))
include(srcdir("single_atom_four_level.jl"))

# load data if `df` not exists
if !isdefined(Main, :df)
    df = collect_results(datadir("zeeman_beat", "gaussian_sampling"))
end

# Make sure if B field noise look like what we expect.
function draw_sampled_Bfield(df::DataFrame)
    # extract Bfields from DataFrame
    Bfields = [df[ii, :params]["Bfield"] for ii in range(1, size(df, 1))]*1e3
    Bx, By, Bz = [[Bfields[ii][jj] for ii in range(1, length(Bfields))] for jj in 1:3]
    # generate scatter plot
    fig1 = scatter(Bx, By, Bz, ms=1, color=:black, lab="")
    fig2 = stephist([norm(Bfields[ii]) for ii in range(1, length(Bfields))], bins=50, normed=true, lab="", xlab="Bfield magnitude [mG]", ylab="Probability")
    fig3 = stephist(Bx, bins=50, normed=true, lab="", xlab="Bx [mG]", ylab="Probability")
    fig4 = stephist(By, bins=50, normed=true, lab="", xlab="By [mG]", ylab="Probability")
    fig5 = stephist(Bz, bins=50, normed=true, lab="", xlab="Bz [mG]", ylab="Probability")
    fig = plot(fig1, fig2, fig3, fig4, fig5, layout=(2, 3), size=(800, 600), suptitle="B-field distribution", titlefontsize=8,)
    return fig
end
fig = draw_sampled_Bfield(df)
savefig(fig, datadir("zeeman_beat", "gaussian_sampling", "Bfield_distribution.pdf"))


function get_avg_fit(df)
    df_result = df[:, :result]
    It_avg = mean([df_result[ii]["It"] for ii in eachindex(df_result)])
    t_out = df_result[1]["t_out"]
    # fit the averaged data
    @. model(x, p) = p[1]*(exp(-(x - p[4])/p[2]) - exp(-(x - p[4])/p[3])) + p[5] 
    p0 = [maximum(It_avg), 10, 1, 0, 0]
    fit_raw = curve_fit(model, t_out, It_avg, p0)
    fit_model = model(t_out, fit_raw.param)
    se = LsqFit.stderror(fit_raw)

    dtauP = (fit_raw.param[2] - 10)/10
    dtauD = (fit_raw.param[3] - 1)

    fig = plot(t_out, It_avg, lab="B-field averaged intensity")
    plot!(t_out, fit_model, lab="double exponential fit", ls=:dash)
    fig_residual = plot(t_out, fit_raw.resid, lab="fit residuals")
    figtot = plot(fig, fig_residual, layout=(2, 1), size=(800, 600), suptitle="ΔτD/τD = $(round(dtauD, digits=4)), ΔτP/τP = $(round(dtauP, digits=3))", titlefontsize=8, xlab="tΓ")
    savefig(figtot, datadir("zeeman_beat", "gaussian_sampling", "avg_fit.pdf"))

    return figtot
end

get_avg_fit(df)