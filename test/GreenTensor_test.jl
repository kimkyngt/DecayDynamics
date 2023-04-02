using DrWatson, Plots
@quickactivate "DecayDynamics"
include("../src/base.jl")


r = range(0, 2, length=1000)
γxx = [ϵ_q(0)*imag(GreenTensor([rr, 0, 0]))*ϵ_q(0)' for rr in r]
γzz = [ϵ_q(0)*imag(GreenTensor([0, 0, rr]))*ϵ_q(0)' for rr in r]
fig = plot(r, real.(γxx), label="γ_xx")
plot!(r, real.(γzz), label="γ_zz", xlab="r", title="k = 2π")
fig
# plot(r, [ϵ_q(0)*real(GreenTensor([0, 0, rr]))*ϵ_q(0)' for rr in r])
# savefig(fig, plotsdir("GreenTensor_test.pdf"))
