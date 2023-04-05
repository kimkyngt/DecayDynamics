using DrWatson, Plots
@quickactivate "DecayDynamics"
include("../src/base.jl")

fig = plot()
r = range(0, 2, length=1000)
γxx = [imag(GreenTensor([rr, 0, 0])[1, 1]) for rr in r]
γzz = [imag(GreenTensor([rr, 0, 0])[3, 3]) for rr in r]
plot!(r, real.(γxx), label="γ_xx", ls=:auto,)
plot!(r, real.(γzz), label="γ_zz", ls=:auto, xlab="r", title="k = 2π")
fig
# plot(r, [ϵ_q(0)*real(GreenTensor([0, 0, rr]))*ϵ_q(0)' for rr in r])
savefig(fig, plotsdir("GreenTensor_test.pdf"))
