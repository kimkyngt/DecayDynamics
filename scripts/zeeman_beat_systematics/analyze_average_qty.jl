using DrWatson, DataFrames, Statistics
include(srcdir("Base.jl"))
include(srcdir("single_atom_four_level.jl"))

# load data if `df` not exists
if !isdefined(Main, :df)
    df = collect_results(datadir("zeeman_beat"))
end

function get_avg_fit(rowindx)
    df_results = df[rowindx, :results]
    It_avg = mean([df_results[ii]["It"] for ii in eachindex(df_results)])
    t_out = df_results[1]["t_out"]

    # fit the averaged data
    @. model(x, p) = p[1]*(exp(-(x - p[4])/p[2]) - exp(-(x - p[4])/p[3])) + p[5] 
    p0 = [maximum(It_avg), 10, 1, 0, 0]
    fit_raw = curve_fit(model, t_out, It_avg, p0)
    fit_model = model(t_out, fit_raw.param)
    se = LsqFit.stderror(fit_raw)
    

    fig = plot(t_out, It_avg, lab="B-field averaged intensity")
    plot!(t_out, fit_model, lab="double exponential fit", ls=:dash)
    fig_residual = plot(t_out, fit_raw.resid, lab="fit residuals")
    figtot = plot(fig, fig_residual, layout=(2, 1), size=(800, 600), suptitle="τD/Γ = $(round(fit_raw.param[3], digits=4))($(round(se[3], digits=4))), τP/Γ = $(round(fit_raw.param[2], digits=3))($(round(se[2], digits=4)))", titlefontsize=8, xlab="tΓ")
    savefig(figtot, datadir("zeeman_beat", df[rowindx, :path]*"avg_fit.pdf"))

    dtauP = (fit_raw.param[2] - 10)/10
    dtauD = (fit_raw.param[3] - 1)
    B = df[rowindx, :params]["B"]
    rho0_indx = df[rowindx, :params]["rho0_indx"]
    thetaDet = df[rowindx, :params]["thetaDet"]

    return dtauD, dtauP, B, rho0_indx, thetaDet
end

df_avg_fit = DataFrame(dtauD=[], dtauP=[], B=[], rho0_indx=[], thetaDet=[])
for rowindx in range(1, size(df, 1))
    dtauD, dtauP, B, rho0_indx, thetaDet = get_avg_fit(rowindx)
    push!(df_avg_fit, [dtauD, dtauP, B, rho0_indx, thetaDet])
end

# find unique value of thetaDet and rho0_indx
thetaDet_list = unique(df_avg_fit[:, "thetaDet"])
rho0_indx_list = unique(df_avg_fit[:, "rho0_indx"])
B_list = unique(df_avg_fit[:, "B"])

function Bfield_dependence(thetDet, rho_indx)
    # fixed theta, rho0_indx and vary B
    _df_avg = filter(row -> row[:thetaDet] == thetDet && row[:rho0_indx] == rho_indx, df_avg_fit)
    fig = plot(_df_avg[:, "B"]*1e3, _df_avg[:, "dtauP"], lab="P",  st=:scatter, ms=5)
    plot!(_df_avg[:, "B"]*1e3, _df_avg[:, "dtauD"], lab="D", ylab="Fractional shift", st=:scatter, title="thetaDet = $(round(thetDet/π*180)) deg, rho0_indx = $(rho_indx)", xlab="B-field magnitude [mG]", ms=5)
    
    savefig(fig, datadir("zeeman_beat", "analyze_average_qty", "Bfield_dependence_thetaDet=$(thetDet)_rho0_indx=$(rho_indx).pdf"))
end

# run the analysis for the B-field dependence
for thetaDet in thetaDet_list
    for rho_indx in rho0_indx_list
        Bfield_dependence(thetaDet, rho_indx)
    end
end


function thetaDet_dependence(rho_indx, B)
    # fixed rho_indx, B and vary thetaDet
    _df_avg = filter(row -> row[:B] == B && row[:rho0_indx] == rho_indx, df_avg_fit)
    fig = plot(_df_avg[:, "thetaDet"]/π*180, _df_avg[:, "dtauP"], lab="P",  st=:scatter, ms=5)
    plot!(_df_avg[:, "thetaDet"]/π*180, _df_avg[:, "dtauD"], lab="D", ylab="Fractional shift", st=:scatter, title="B = $(round(B*1e3)) mG, rho0_indx = $(rho_indx)", xlab="thetaDet [deg]", ms=5)
    
    savefig(fig, datadir("zeeman_beat", "analyze_average_qty", "thetaDet_dependence_B=$(round(B*1e3))_rho0_indx=$(rho_indx).pdf"))
end

# run the analysis for the thetaDet dependence
for B in B_list
    for rho_indx in rho0_indx_list
        thetaDet_dependence(rho_indx, B)
    end
end