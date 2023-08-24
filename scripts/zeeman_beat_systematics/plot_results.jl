using DrWatson, DataFrames
include(srcdir("Base.jl"))
include(srcdir("single_atom_four_level.jl"))

# load data if `df` not exists
if !isdefined(Main, :df)
    df = collect_results(datadir("zeeman_beat"))
end

function plot_results(df, rownum)
    results = df[rownum, :results]
    filename = df[rownum, :path]
    # Plot the results
    thetas = []
    phis = []
    tauPs = []
    tauDs = []
    for ii in eachindex(results)
        Bfield = results[ii]["Bfield"]
        fit_raw = results[ii]["fit_raw"]

        theta, phi = xyz_2_polar(Bfield[1], Bfield[2], Bfield[3])
        push!(thetas, theta)
        push!(phis, phi)
        push!(tauPs, fit_raw.param[2])
        push!(tauDs, fit_raw.param[3])
    end

    dtauP = (tauPs .- 10)/10
    dtauPmax = maximum(abs.(dtauP))
    dtauD = (tauDs .- 1)
    dtauDmax = maximum(abs.(dtauD))

    fig = plot(
        scatter(thetas/π, phis/π, marker_z=dtauP, markersize=11, markerstrokewidth=0, lab="τP fractional error", xlab="θ/π", ylab="ϕ/π", legend=false, clim=(-dtauPmax, dtauPmax), colorbar=false, color=:vik, title="P"),
        scatter(thetas/π, phis/π, marker_z=dtauD, markersize=11, markerstrokewidth=0, lab="τD, fractional error", xlab="θ/π", ylab="ϕ/π", legend=false, clim=(-dtauDmax, dtauDmax), colorbar=false, color=:vik, title="D"),
        plot(phis/π, dtauP, marker_z=dtauP, st=:scatter, markersize=11, markerstrokewidth=0, ylab="fractional error", xlab="ϕ/π", legend=false, clim=(-dtauPmax, dtauPmax), color=:vik, cb=:none, ylim=[-dtauPmax, dtauPmax]*1.1),
        plot(phis/π, dtauD, marker_z=dtauD, st=:scatter, markersize=11, markerstrokewidth=0, ylab="fractional error", xlab="ϕ/π", legend=false, clim=(-dtauDmax, dtauDmax), color=:vik, cb=:none, ylim=[-dtauDmax, dtauDmax]*1.1),
    size=(800, 800), suptitle=splitpath(filename)[end], )
    return fig
end

# save the figures for each row of `df`
for ii in range(1, size(df, 1))
    fig = plot_results(df, ii)
    savefig(fig, datadir("zeeman_beat", df[ii, :path]*"fig.pdf"))
end


