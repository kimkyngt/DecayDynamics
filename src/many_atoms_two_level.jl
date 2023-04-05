using Plots, QuantumOptics, LaTeXStrings, WignerSymbols, ProgressMeter, Latexify

# """
# Σ_iq for two atoms.
# """
# function Σ_iq(i::Int, q::Int, F_1, F_2)
#     comp_basis = SpinBasis(F_i[1]) ⊕ SpinBasis(F_i[2])    
#     if i == 1
#         return Σ_q(q, F_1, F_2) ⊗ identityoperator(comp_basis)    
#     elseif i == 2
#         return identityoperator(comp_basis) ⊗ Σ_q(q, F_1, F_2)
#     else
#         throw(ArgumentError("Argument i must be one of [1, 2]"))
#     end
# end


"""
Solve time dynamics with master equation.
"""
function evolve_master(parameters::Dict; check_rates::Bool=false)
    @unpack F_i, g_i, Γ_i, Bfield, positions, ρ_0, tspan = parameters

    J(i, j, p, q) = - 3/4*Γ_i[1]*ϵ_q(p)*real.(GreenTensor(positions[i] - positions[j])) * adjoint(ϵ_q(q))
    Γ(i, j, p, q) = 3/2*Γ_i[1]*ϵ_q(p)*imag.(GreenTensor(positions[i] - positions[j])) * adjoint(ϵ_q(q))

    q_list = [-1, 0, 1]
    ptlindx = range(1, length(positions))
    
    H = sum(J(i, j, p, q)  * (dagger(Σ_iq(i, p, F_i))*Σ_iq(j, q, F_i)) for i in ptlindx for j in ptlindx for p in q_list for q in q_list)

    indx = [(i, p) for i in ptlindx for p in q_list]
    # Jump = [Σ_iq(i, p, F_i) for (i, p) in indx]
    # rates = [Γ(i, j, p, q)  for (i, p) in indx, (j, q) in indx]
    Jump        = [Σ_iq(i, p, F_i)           for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    Jump_dagger = [dagger(Σ_iq(j, q, F_i))  for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    rates       = [Γ(i, j, p, q)                  for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    Jmatrix     = [J(i, j, p, q)                  for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]

    if check_rates
        # println("Re[Gamma]:")
        # display(real.(rates))
        # println("Basis, (particle, polarization):")
        # display(indx)
        println("Re[Omega]:")
        display(real.(Jmatrix))
    end

    p = Progress(length(tspan), "Solving master equation...")
    function forward_progress(t, rho) 
        # function for checking the progress
        next!(p)
        return copy(rho)
    end
    t_out, ρ_t = timeevolution.master(tspan, ρ_0, H, Jump, 
        rates=rates, 
        Jdagger=Jump_dagger,
        fout=forward_progress,
    )
    result = copy(parameters)
    @pack! result = t_out, ρ_t, Jump, rates
    return result
end

"""
Plot time evolution data
"""
function plot_dynamics(result::Dict; kwargs...)
    @unpack ρ_t, t_out, F_i, positions, Γ_i, m_exc= result
    # Plot population dynamics
    rhoee_1 = [real.(expect( (identityoperator(SpinBasis(F_i[1])) ⊕ 0*identityoperator(SpinBasis(F_i[2]))), ptrace(ρ_t[ii], 2) ) ) for ii in eachindex(ρ_t)]
    rhogg_1 = [real.(expect( (0*identityoperator(SpinBasis(F_i[1])) ⊕ identityoperator(SpinBasis(F_i[2]))), ptrace(ρ_t[ii], 2) ) ) for ii in eachindex(ρ_t)]
    rhoee_2 = [real.(expect( (identityoperator(SpinBasis(F_i[1])) ⊕ 0*identityoperator(SpinBasis(F_i[2]))), ptrace(ρ_t[ii], 1) ) ) for ii in eachindex(ρ_t)]
    rhogg_2 = [real.(expect( (0*identityoperator(SpinBasis(F_i[1])) ⊕ identityoperator(SpinBasis(F_i[2]))), ptrace(ρ_t[ii], 1) ) ) for ii in eachindex(ρ_t)]
    fig1 = plot(title="Spin: $(latexify(F_i[1])) → $(latexify(F_i[2])), "*L"\vec{r}_{12} = "*"$(round.(positions[1]-positions[2], digits=2))", xlab="tΓ", leg=:outerright)
    plot!(fig1, t_out*Γ_i[1], rhoee_1, label="ptl1, e", ls=:solid, lw=2)
    plot!(fig1, t_out*Γ_i[1], rhogg_1, label="ptl1, g", ls=:solid, lw=2)
    plot!(fig1, t_out*Γ_i[1], rhoee_2, label="ptl2, e", ls=:solid, lw=2)
    plot!(fig1, t_out*Γ_i[1], rhogg_2, label="ptl2, g", ls=:solid, lw=2)
    plot!(fig1, t_out*Γ_i[1], real.(tr.(ρ_t)), label="Tr(ρ)", ls=:dash, lc=:black)
    plot!(tspan*Γ_i[1], exp.(-Γ_i[1]*tspan), lab=L"\propto e^{-\Gamma t}", color=:black, ls=:dash, lw=1)
    plot!(tspan*Γ_i[1], exp.(-(2)*Γ_i[1]*tspan), lab=L"\propto e^{-2\Gamma t}", color=:blue, ls=:dash, lw=1)
    plot!(tspan*Γ_i[1], 1 .- exp.(-Γ_i[1]*tspan), lab=L"\propto 1-e^{-\Gamma t}", color=:black, ls=:dash, lw=1)
    plot!(tspan*Γ_i[1], 1 .- exp.(-(2)*Γ_i[1]*tspan), lab=L"\propto 1-e^{-2\Gamma t}", color=:red, ls=:dash, lw=1)
    # plot!(leg=:outerbottom)

    # two particle populations
    rhoee =  [real.(expect( projector(normalize( (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2]))) ⊗ ((Fm_state(F_i[1], m_exc)) ⊕ Ket(SpinBasis(F_i[2])) ))), ρ_t[ii] ) ) for ii in eachindex(ρ_t)]
    rhogg =  [real.(expect( projector(normalize( (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], 0)) ⊗ (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], 0)) )),   ρ_t[ii] ) ) for ii in eachindex(ρ_t)]
    rhoss =  [real.(expect( projector(normalize( (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2]))) ⊗ (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], 0)) + 
        (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], 0)) ⊗ (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2]))) )), ρ_t[ii] ) ) for ii in eachindex(ρ_t)]
    rhoas =  [real.(expect( projector(normalize( (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2]))) ⊗ (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], 0)) - 
        (Ket(SpinBasis(F_i[1])) ⊕ Fm_state(F_i[2], 0)) ⊗ (Fm_state(F_i[1], m_exc) ⊕ Ket(SpinBasis(F_i[2]))) )), ρ_t[ii] ) ) for ii in eachindex(ρ_t)]
    C = 2*[abs.(expect(
        σ(F_i[1], m_exc, F_i[2], 0) ⊗ dagger(σ(F_i[1], m_exc, F_i[2], 0))
        , ρ_t[ii] ) ) for ii in eachindex(ρ_t)] # Concurrence
    fig2 = plot(xlab="tΓ", leg=:outerright)
    plot!(fig2, t_out*Γ_i[1], rhoee, label="ee",)
    plot!(fig2, t_out*Γ_i[1], rhogg, label="gg",)
    plot!(fig2, t_out*Γ_i[1], rhoss, label="eg+ge",)
    plot!(fig2, t_out*Γ_i[1], rhoas, label="eg-ge",)
    plot!(fig2, t_out*Γ_i[1], C    , label="Concurrence")
    plot!(tspan*Γ_i[1], rhoee[1]*exp.(-Γ_i[1]*tspan), lab=L"\propto e^{-\Gamma t}", color=:black, ls=:dash, lw=1)
    plot!(tspan*Γ_i[1], rhoee[1]*exp.(-(2)*Γ_i[1]*tspan), lab=L"\propto e^{-2\Gamma t}", color=:blue, ls=:dash, lw=1)
    plot!(tspan*Γ_i[1], 1 .- exp.(-Γ_i[1]*tspan), lab=L"\propto 1-e^{-\Gamma t}", color=:black, ls=:dash, lw=1)
    plot!(tspan*Γ_i[1], 1 .- exp.(-(2)*Γ_i[1]*tspan), lab=L"\propto 1-e^{-2\Gamma t}", color=:red, ls=:dash, lw=1, yticks=0:0.2:1)

    # Position of the atoms
    fig3 = plot(xlabel="x", ylabel="y", zlabel="z", title="Atom positions", aspect_ratio=:equal)
    for ii in eachindex(positions)
        plot!(fig3, [positions[ii][1]], [positions[ii][2]], [positions[ii][3]], st=:scatter, label="Atom $(ii)")
    end

    fig = plot(fig1, fig2, layout=(2, 1), size=(800, 600))
    return fig
end