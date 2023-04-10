using Distributions, Plots

μ = [0.0, 0.0, 0.0]
σ = [1, 1, 100]
Σ = [σ[1]^2 0.0 0.0; 0.0 σ[2]^2 0.0; 0.0 0.0 σ[3]^2]
dist = MvNormal(μ, Σ)
distance_dist = MvNormal(μ, 2*Σ)
Nsample = 10_000
x_pdf = range(-1, 1, length=100)

r = rand(dist, Nsample)
fig1 = []
for ii in [1, 2, 3]
    global fig1
    fig = plot(r[ii, :], st=:histogram, legend=false, normalize=:pdf)
    xx = x_pdf*5*σ[ii]
    plot!(fig, xx, 2π*σ[mod1(ii+1, 3)]*σ[mod1(ii+2, 3)]*[pdf(dist, [jj == ii ? x : 0 for jj in [1, 2, 3]]) for x in xx ], label="pdf",  lw=2)
    push!(fig1, fig)
end
fig1 = plot(fig1..., plot_title = "Atom sampling", layout=(3, 1), size= (400, 600))

d = Matrix(undef,3, Nsample)
d = rand(dist, Nsample) - rand(dist, Nsample)
fig2 = []
for ii in [1, 2, 3]
    global fig1
    fig = plot(d[ii, :], st=:histogram, legend=false, normalize=:pdf)
    xx = x_pdf*5*σ[ii]*sqrt(2)
    plot!(fig, xx, 2π*2*σ[mod1(ii+1, 3)]*σ[mod1(ii+2, 3)]*[pdf(distance_dist, [jj == ii ? x : 0 for jj in [1, 2, 3]]) for x in xx ], label="pdf",  lw=2)
    push!(fig2, fig)
end
fig2 = plot(fig2..., plot_title = "Distance sampling", layout=(3, 1), size= (400, 600))

plot(fig1, fig2, layout=(1, 2), size= (800, 600))

