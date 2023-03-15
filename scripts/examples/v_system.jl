using QuantumOptics
using Plots, LaTeXStrings

# Energy level structure
# 2,3 -- > 1

# Parameters
γ_21 = 1
γ_31 = 1

ω3 = 2
ω2 = 0

tmax = 5
dt = 0.002
tspan = 0:dt:tmax

# Basis and operators
b = NLevelBasis(3)

J_12 = transition(b, 1, 2)
J_13 = transition(b, 1, 3)
J_23 = transition(b, 2, 3)

proj_1 = transition(b, 1, 1)
proj_2 = transition(b, 2, 2)
proj_3 = transition(b, 3, 3)

H       = dense(ω3*proj_3 + ω2*proj_2)
J       = [J_12, J_13]
rates   = [γ_21, γ_31]
j       = (sqrt.(rates) .* J)
jdagger = dagger.(j)
J_final = [sum(j[1:2])] 
# J_final = [j[1], j[2]]

# Solve
ψ_0 = 0*nlevelstate(b, 2) + nlevelstate(b, 3)
normalize!(ψ_0)
t_out, rho_t = timeevolution.master_h(
    tspan, ψ_0, H, J_final, 
    # jdagger=J_final_dagger
    )

fig1 = plot()
plot!(fig1, t_out, [real(expect(proj_3, rho_t[ii])) for ii ∈ range(1, length(rho_t))], label=L"\rho^{33}")
plot!(fig1, t_out, [real(expect(proj_2, rho_t[ii])) for ii ∈ range(1, length(rho_t))], label=L"\rho^{22}")
# plot!(fig1, t_out, [real(expect(proj_2 + proj_3, rho_t[ii])) for ii ∈ range(1, length(rho_t))], label=L"\rho^{33} + \rho^{22}")
plot!(fig1, t_out, [real(expect(proj_1, rho_t[ii])) for ii ∈ range(1, length(rho_t))], label=L"\rho^{11}")
plot!(fig1, t_out, [real(expect(one(b), rho_t[ii])) for ii ∈ range(1, length(rho_t))], label=L"Tr(ρ)", lc=:black)

fig2 = plot()
plot!(fig2, t_out, [real(expect(J_23, rho_t[ii])) for ii ∈ range(1, length(rho_t))], label=L"\rho^{23}")
plot!(fig2, t_out, [real(expect(J_12, rho_t[ii])) for ii ∈ range(1, length(rho_t))], label=L"\rho^{12}")

plot(
    fig1, fig2, size=(600, 600), layout=(2, 1), 
    legend=:right,
    xlab="time"
    )
