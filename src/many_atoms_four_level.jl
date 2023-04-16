using Plots, QuantumOptics, LaTeXStrings, WignerSymbols, ProgressMeter, Latexify, LsqFit

"""
Solve dynamics with master equation.
"""
function evolve_master(parameters::Dict; check_rates::Bool=false)
    @unpack F_i, Γ_kl, positions, ρ_0, tspan, knorm,  = parameters

    # J and Γ matrix
    J(i, j, p, q, kl::CartesianIndex) = - 3/4*Γ_kl[kl]*ϵ_q(p)*real(GreenTensor(positions[i] - positions[j], knorm[kl])) * adjoint(ϵ_q(q))
    Γ(i, j, p, q, kl::CartesianIndex) =   3/2*Γ_kl[kl]*ϵ_q(p)*imag(GreenTensor(positions[i] - positions[j], knorm[kl])) * adjoint(ϵ_q(q))

    q_list = [-1, 0, 1]
    ptlindx = range(1, length(positions))
    transition_indx = findall(!iszero, Γ_kl)

    H = sparse(sum(J(i, j, p, q, kl) * (dagger(Σ_iq(i, p, F_i, kl[1], kl[2]))*Σ_iq(j, q, F_i, kl[1], kl[2])) for i in ptlindx for j in ptlindx for p in q_list for q in q_list for kl in transition_indx))

    Jump = [Σ_iq(i, p, F_i, kl[1], kl[2])                           for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list for kl in transition_indx]
    Jump_dagger = [dagger(Σ_iq(j, q, F_i, kl[1], kl[2]))    for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list for kl in transition_indx]
    rates = [Γ(i, j, p, q, kl)                                      for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list for kl in transition_indx]
 
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
        abstol=1e-9,
        reltol=1e-9,
    )
    result = copy(parameters)
    @pack! result = t_out, ρ_t
    return result
end

"""
Plot time evolution data
"""
function plot_dynamics(result::Dict; return_dict::Bool=false, kwargs...)
    @unpack ρ_t, t_out, F_i, positions, Γ_kl, m_exc, exc_frac = result
    
    # state_label = ["³D₁", "³P₀", "³P₁" "¹S₀"]
    state_label = ["e", "f1", "f0", "g"]
    τtot = 1/(Γ_kl[1, 3] + Γ_kl[1, 2])
    τ1 = 1/Γ_kl[2, 4]
    τ2 = 1/Γ_kl[1, 2]
    tΓ = t_out/τtot#*Γ_kl[1, 2]
    # Plot population dynamics

    fig1 = plot(xlab="tΓ", leg=:right)
    for ii in eachindex(state_label)
        # plot!(
        #     fig1, 
        #     tΓ,
        #     [real.(expect(⊕([kk == ii ? identityoperator(SpinBasis(F_i[kk])) : projector(Ket(SpinBasis(F_i[kk]))) for kk in eachindex(F_i)]...), ptrace(ρ_t[jj], 2)) ) for jj in eachindex(ρ_t)],
        #     label="ptl1, $(state_label[ii])",
        # )
        plot!(
            fig1, 
            tΓ,
            [real.(expect(  
                ⊕([kk == ii ? identityoperator(SpinBasis(F_i[kk])) : projector(Ket(SpinBasis(F_i[kk]))) for kk in eachindex(F_i)]...) 
                ⊗ 
                ⊕([identityoperator(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
                , ρ_t[jj]) ) for jj in eachindex(ρ_t)],
            label="$(state_label[ii]) population",
            ls=:solid
        )
        # plot!(
        #     fig1, 
        #     tΓ,
        #     [real.(expect( ⊕( (circshift([ii == 1 ? 1 : 0 for ii in eachindex(F_i)], ii-1) .* [identityoperator(SpinBasis(F)) for F in F_i])...), ptrace(ρ_t[jj], 1)) ) for jj in eachindex(ρ_t)],
        #     label="ptl2, $(state_label[ii])",
        # )
    end
    plot!(fig1, tΓ, real.(tr.(ρ_t)), label="Tr(ρ)", ls=:dash, lc=:black)
    analytic_model = τtot*τ1/((τ1-τtot)*τ2)*(exp.(-t_out/τ1) .- exp.(-t_out/τtot))
    # plot!(fig1, tΓ, analytic_model, lab="Analytic", color=:black, ls=:dash)
    # plot!(fig1, tΓ, exp.(-t_out/τtot), lab=L"\propto e^{-\Gamma_{e\rightarrow f} t}", color=:blue, ls=:dash)
    # plot!(fig1, tΓ, exp.(-2t_out/τtot), lab=L"\propto e^{-2\Gamma_{e\rightarrow f} t}", color=:red, ls=:dash)

    # two particle populations
    fig2 = plot(xlab="tΓ", leg=:right)
    for ii in eachindex(state_label)
        # # Diagonal 
        # plot!(
        #     fig2, 
        #     tΓ,
        #     [real.(expect( projector( normalize(
        #         ⊕( (circshift([ii == 1 ? 1 : 0 for ii in eachindex(F_i)], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...) ⊗ ⊕( (circshift([ii == 1 ? 1 : 0 for ii in eachindex(F_i)], ii-1) .* [spinup(SpinBasis(F)) for F in F_i])...)
        #         )), ρ_t[jj]) ) for jj in eachindex(ρ_t)],
        #     label="$(state_label[ii])$(state_label[ii])",
        # )
        
        # Symmetric coherence
        plot!(
            fig2, 
            tΓ,
            [real.(expect( projector( normalize(
                ⊕([kk == ii ? spinup(SpinBasis(F_i[kk])) : Ket(SpinBasis(F_i[kk])) for kk in eachindex(F_i)]...) ⊗ ⊕( [kk == mod1(ii+1, length(F_i)) ? spinup(SpinBasis(F_i[kk])) : Ket(SpinBasis(F_i[kk])) for kk in eachindex(F_i)]...) + 
                ⊕( [kk == mod1(ii+1, length(F_i)) ? spinup(SpinBasis(F_i[kk])) : Ket(SpinBasis(F_i[kk])) for kk in eachindex(F_i)]...) ⊗ ⊕( ([kk == ii ? spinup(SpinBasis(F_i[kk])) : Ket(SpinBasis(F_i[kk])) for kk in eachindex(F_i)]...)))), ρ_t[jj]) ) for jj in eachindex(ρ_t)],
            label="$(state_label[ii])$(state_label[mod1(ii+1, length(F_i))]) + $(state_label[mod1(ii+1, length(F_i))])$(state_label[ii])",
        )

        # Antisymmetric coherence
        plot!(
            fig2, 
            tΓ,
            [real.(expect( projector( normalize(
                ⊕([kk == ii ? spinup(SpinBasis(F_i[kk])) : Ket(SpinBasis(F_i[kk])) for kk in eachindex(F_i)]...) ⊗ ⊕( [kk == mod1(ii+1, length(F_i)) ? spinup(SpinBasis(F_i[kk])) : Ket(SpinBasis(F_i[kk])) for kk in eachindex(F_i)]...) - 
                ⊕( [kk == mod1(ii+1, length(F_i)) ? spinup(SpinBasis(F_i[kk])) : Ket(SpinBasis(F_i[kk])) for kk in eachindex(F_i)]...) ⊗ ⊕( ([kk == ii ? spinup(SpinBasis(F_i[kk])) : Ket(SpinBasis(F_i[kk])) for kk in eachindex(F_i)]...)))), ρ_t[jj]) ) for jj in eachindex(ρ_t)],
            label="$(state_label[ii])$(state_label[mod1(ii+1, length(F_i))]) - $(state_label[mod1(ii+1, length(F_i))])$(state_label[ii])",
            ls=:dash
        )
  
    end
    # Comparison to analytic and fitting
    xdata = t_out
    ydata = [real.(expect(⊕([kk == 2 ? identityoperator(SpinBasis(F_i[kk])) : projector(Ket(SpinBasis(F_i[kk]))) for kk in eachindex(F_i)]...), ptrace(ρ_t[jj], 2)) ) for jj in eachindex(ρ_t)]

    # fit data
    @. model(x, p) = p[1]*(exp(-(x - p[4])/p[2]) - exp(-(x - p[4])/p[3])) + p[5]
    transition_indx = findall(!iszero, Γ_kl)
    p0 = [maximum(ydata), 1/Γ_kl[transition_indx[3]], τtot, 0, 0] # True taus
    fit_raw = curve_fit(model, xdata, ydata, p0)
    fit_model = model(xdata, fit_raw.param)

    fig3 = plot(xlab="tΓ",)
    plot!(fig3, tΓ, ydata, label="data",)
    plot!(fig3, tΓ, fit_model, label="fit",)

    fig4 = plot(xlab="tΓ",)
    plot!(fig4, tΓ, (ydata - fit_model)/maximum(fit_model), lab="(data - fit)/max[fit], \n (τᴾ/τᴾtrue, τᴰ/τᴰtrue) = ($(round(fit_raw.param[2]/p0[2], digits=4)), $(round(fit_raw.param[3]/p0[3], digits=4)))", )

    fig = plot(fig1, fig3, fig2, fig4, layout=(2, 2), size=(1200, 800), plot_title=" $(state_label[1])($(latexify(F_i[1])), $(latexify(m_exc))) → $(state_label[2])($(latexify(F_i[2]))), $(state_label[3])($(latexify(F_i[3]))) → $(state_label[4])($(latexify(F_i[4]))), "*L"\vec{r}_{12} = "*"$(round.(positions[1]-positions[2], digits=2)). Exc. frac: $(round(exc_frac[1], digits=3))", margins=5Plots.mm)

    # Saving data
    chisq = sum((ydata - fit_model).^2)
    tauDratio = fit_raw.param[3]/p0[3]
    tauPratio = fit_raw.param[2]/p0[2]

    fitresult = @strdict chisq tauDratio tauPratio fit_raw

    if return_dict
        return fig, fitresult
    else
        return fig
    end
end