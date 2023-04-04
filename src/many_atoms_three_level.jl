using Plots, QuantumOptics, LaTeXStrings, WignerSymbols, ProgressMeter, Latexify

"""
Solve time dynamics with master equation.
"""
function evolve_master(parameters::Dict; check_rates::Bool=false)
    @unpack F_i, g_i, Γ_i, Bfield, positions, ρ_0, tspan, knorm = parameters

    # J and Γ mattrix
    J(i, j, p, q, k) = - 3/4*Γ_i[k]*ϵ_q(p)*real.(GreenTensor(positions[i] - positions[j], knorm[k])) * adjoint(ϵ_q(q))
    Γ(i, j, p, q, k) =   3/2*Γ_i[k]*ϵ_q(p)*imag.(GreenTensor(positions[i] - positions[j], knorm[k])) * adjoint(ϵ_q(q))

    q_list = [-1, 0, 1]
    ptlindx = range(1, length(positions))
    
    H = sum(J(i, j, p, q, k)  * (dagger(Σ_iq(i, p, F_i, k, k+1))*Σ_iq(j, q, F_i, k, k+1)) for i in ptlindx for j in ptlindx for p in q_list for q in q_list for k in 1:(length(knorm)-1))

    indx = [(i, p) for i in ptlindx for p in q_list]
    Jump        = [Σ_iq(i, p, F_i, k, k+1)               for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list for k in 1:(length(knorm)-1)]
    Jump_dagger = [dagger(Σ_iq(j, q, F_i, k, k+1))       for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list for k in 1:(length(knorm)-1)]
    rates       = [Γ(i, j, p, q, k)                       for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list for k in 1:(length(knorm)-1)]
 
    if check_rates
        println("Re[Gamma]:")
        display(real.(rates))
        println("Basis, (particle, polarization):")
        display(indx)
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
    @pack! result = t_out, ρ_t, rates
    return result
end

"""
Plot time evolution data
"""
function plot_dynamics(result::Dict; kwargs...)
    @unpack ρ_t, t_out, F_i, positions, Γ_i, m_exc= result
    
    state_label = ["e", "f", "g"]
    tΓ = t_out*Γ_i[1]
    # Plot population dynamics

    fig1 = plot(title="Spin: ($(latexify(F_i[1])), $(latexify(m_exc))) → $(latexify(F_i[2])) → $(latexify(F_i[3])), "*L"\vec{r}_{12} = "*"$(round.(positions[1]-positions[2], digits=2))", xlab="tΓ₁", leg=:outerright)
    for ii in eachindex(state_label)
        plot!(
            fig1, 
            tΓ,
            [real.(expect( ⊕( (circshift([1, 0, 0], ii-1) .* [identityoperator(SpinBasis(F)) for F in F_i])...), ptrace(ρ_t[jj], 2)) ) for jj in eachindex(ρ_t)],
            label="ptl1, $(state_label[ii])",
        )
        plot!(
            fig1, 
            tΓ,
            [real.(expect( ⊕( (circshift([1, 0, 0], ii-1) .* [identityoperator(SpinBasis(F)) for F in F_i])...), ptrace(ρ_t[jj], 1)) ) for jj in eachindex(ρ_t)],
            label="ptl2, $(state_label[ii])",
        )
    end
    plot!(fig1, tΓ, real.(tr.(ρ_t)), label="Tr(ρ)", ls=:dash, lc=:black)
    plot!(fig1, tΓ, Γ_i[1]/(Γ_i[1] - Γ_i[2])*(exp.(-Γ_i[2]*t_out) .- exp.(-Γ_i[1]*t_out)), lab="Analytic", color=:black, ls=:dash, lw=1)
    plot!(fig1, tΓ, exp.(-Γ_i[1]*t_out), lab=L"\propto e^{-\Gamma_{D} t}", color=:blue, ls=:dash, lw=1)
    plot!(fig1, tΓ, exp.(-2Γ_i[1]*t_out), lab=L"\propto e^{-2\Gamma_{D} t}", color=:red, ls=:dash, lw=1)

    # two particle populations
    fig2 = plot(xlab="tΓ₁", leg=:outerright)
    for ii in eachindex(state_label)
        plot!(
            fig2, 
            tΓ,
            [real.(expect( projector( normalize(
                ⊕( (circshift([1, 0, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...) ⊗ ⊕( (circshift([1, 0, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...)
                )), ρ_t[jj]) ) for jj in eachindex(ρ_t)],
            label="$(state_label[ii])$(state_label[ii])",
        )
        plot!(
            fig2, 
            tΓ,
            [real.(expect( projector( normalize(
                ⊕( (circshift([1, 0, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...) ⊗ ⊕( (circshift([0, 1, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...) + 
                ⊕( (circshift([0, 1, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...) ⊗ ⊕( (circshift([1, 0, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...)
                )), ρ_t[jj]) ) for jj in eachindex(ρ_t)],
            label="$(state_label[ii])$(state_label[mod1(ii+1, 3)]) + $(state_label[mod1(ii+1, 3)])$(state_label[ii])",
        )
        plot!(
            fig2, 
            tΓ,
            [real.(expect( projector( normalize(
                ⊕( (circshift([1, 0, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...) ⊗ ⊕( (circshift([0, 1, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...) - 
                ⊕( (circshift([0, 1, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...) ⊗ ⊕( (circshift([1, 0, 0], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...)
                )), ρ_t[jj]) ) for jj in eachindex(ρ_t)],
            label="$(state_label[ii])$(state_label[mod1(ii+1, 3)]) - $(state_label[mod1(ii+1, 3)])$(state_label[ii])",
            ls=:dash
        )
    end

    fig = plot(fig1, fig2, layout=(2, 1), size=(800, 600))
    return fig
end