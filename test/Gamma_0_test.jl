using DrWatson, Plots
@quickactivate "DecayDynamics"
include("../src/base.jl")

Fg = 5//1
Fe = 6//1

C(m, q) = (-1)^(Fg - m) * wigner3j(Fg, 1, Fe, -m, q, m-q)
Gam(Fe, Fg) = (2*Fe + 1)


for me = -Fe:Fe
    println("me = $me")
    println("\n Fe, Fg = ($(Fe), $(Fg))")
    println("Fe dependent sum:")
    println(sum(C(me+q, q)^2 for q = [-1, 0, 1] if abs(me + q) <= Fg))
    println("Fe independent sum:")
    println(sum(Gam(Fe, Fg)*C(me+q, q)^2 for q = [-1, 0, 1] if abs(me + q) <= Fg))
end
