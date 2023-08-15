using Plots, QuantumOptics, LaTeXStrings, WignerSymbols, ProgressMeter, Latexify, LsqFit

"""
Solve dynamics with master equation.
parameters:
    - F_i: array of F quantum numbers 
    - g_i: array of g factors
    - Γ_kl: matrix of decay rates with indices (k, l) corresponding to (F_i[k], F_i[l])
    - Bfield: magnetic field strength in units of 
    - ρ_0: initial density matrix
    - tspan: time span to output
"""
function evolve_master(parameters::Dict;)
    @unpack F_i, g_i, Γ_kl, Bfield, ρ_0, tspan = parameters
    # Define operators
    q_list = [-1, 0, 1]
    transition_indx = findall(!iszero, Γ_kl)

    # Construct Hamiltonian
    ## Zeeman shift
    H  = directsum([Bfield[1]*g_i[ii]*F_i[ii]*sigmax(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
    H += directsum([Bfield[2]*g_i[ii]*F_i[ii]*sigmay(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
    H += directsum([Bfield[3]*g_i[ii]*F_i[ii]*sigmaz(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
    H = 2π*sparse(H)

    # Jump operators
    Jump = [Σ_q(q, F_i, kl[1], kl[2]) for q in q_list for kl in transition_indx]
    rates = [Γ_kl[kl[1], kl[2]] for q in q_list for kl in transition_indx]
 
    p = Progress(length(tspan), "Solving master equation...")
    function forward_progress(t, rho) 
        # function for checking the progress
        next!(p)
        return copy(rho)
    end
    t_out, ρ_t = timeevolution.master(tspan, ρ_0, H, Jump, 
        rates=rates, 
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
function plot_dynamics(result::Dict; kwargs...)
    @unpack F_i, g_i, Γ_kl, Bfield, ρ_t, t_out = result

    b_i = [SpinBasis(F_i[ii]) for ii in eachindex(F_i)]
    
    figs = []
    # Plot population
    for indx in eachindex(F_i)
        fig = plot(ylab="Population, Level-$indx, $(F_i[indx])", legend=:best, xlab="tΓ")
        id_ops = [one(b_i[ii]) for ii in eachindex(F_i)]
        id_vec = zeros(length(F_i)); id_vec[1] = 1; id_vec = circshift(id_vec, indx-1)

        # Plot population of level-indx
        plot!(fig, t_out, real.(expect(directsum([one(b_i[ii]) for ii in eachindex(F_i)]...), ρ_t)), label="Tr(ρ)", ls=:solid, lc=:black)
        plot!(fig, t_out, real.(expect(directsum([id_vec[ii]*id_ops[ii] for ii in eachindex(F_i)]...), ρ_t)), label="Tr(ρ$indx)",  ls=:dash, lc=:black)

        # Plot each sublevel's population
        mycolor = cgrad(:Spectral, Int(2*F_i[indx]+1), categorical=true)
        blank_state = [Ket(b_i[ii]) for ii in eachindex(F_i)]
        for ii in range(1, Int(F_i[indx]*2+1))
            spin = F_i[indx] - (ii-1)
            _blank_state = copy(blank_state)
            _blank_state[indx] = Fm_state(F_i[indx], spin)
            state_to_plot = directsum(_blank_state...)
            plot!(fig, t_out, real.(expect(projector(state_to_plot), ρ_t)), label="m = $spin", color=mycolor[ii])
        end
        push!(figs, fig)
    end
    for indx in eachindex(F_i)
        fig = plot(ylab="Coherence, Level-$indx, $(F_i[indx])", legend=:best, xlab="tΓ")
        # Plot each sublevel's coherence
        # mycolor = cgrad(:Spectral, Int(2*F_i[indx]+1), categorical=true)
        blank_state = [projector(Ket(b_i[ii])) for ii in eachindex(F_i)]
        for ii in range(1, min(3, Int(F_i[indx]*2)))
            spin_h = F_i[indx] - (ii-1)
            spin_l = F_i[indx] - ii
            _blank_state = copy(blank_state)
            _blank_state[indx] = projector(Fm_state(F_i[indx], spin_h), dagger(Fm_state(F_i[indx], spin_l)))
            state_to_plot = directsum(_blank_state...)
            plot!(fig, t_out, real.(expect(state_to_plot, ρ_t)), label="($spin_h, $spin_l)", 
            # color=mycolor[ii], 
            ylim=(-0.5, 0.5),)
        end
        push!(figs, fig)
    end

    # return plot(figs..., size=(1500, 800), layout=(2, 3), margins=30Plots.px;kwargs...)
    return figs
end