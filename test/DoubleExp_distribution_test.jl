using Plots
include("../src/DoubleExp.jl")

tauP = 21.0
tauD = 2.1

fig = stephist([rand(DoubleExp(tauP, tauD)) for _ in 1:100_000], label="sample", normalize=:pdf, )
t = range(0, 100, length=1000)
plot!(fig, t, pdf.(DoubleExp(tauP, tauD), t), label="pdf", xlim=(-1, 50))
fig 