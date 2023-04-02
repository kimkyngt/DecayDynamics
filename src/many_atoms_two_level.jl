using Plots, QuantumOptics, LaTeXStrings, WignerSymbols, ProgressMeter, Latexify

"""
Σ_iq for two atoms.
"""
function Σ_iq(i::Int, q::Int, F_1, F_2)
    comp_basis = SpinBasis(F_i[1]) ⊕ SpinBasis(F_i[2])    
    if i == 1
        return Σ_q(q, F_1, F_2) ⊗ identityoperator(comp_basis)    
    elseif i == 2
        return identityoperator(comp_basis) ⊗ Σ_q(q, F_1, F_2)
    else
        throw(ArgumentError("Argument i must be one of [1, 2]"))
    end
end


"""
Solve time dynamics with master equation.
"""
function evolve_master(parameters::Dict; check_rates::Bool=false)
    @unpack F_i, g_i, Γ_i, Bfield, positions, ρ_0, tspan = parameters

    J(i, j, p, q) = - 3/4*Γ_i[1]*ϵ_q(p)*real.(GreenTensor(positions[i] - positions[j])) * adjoint(ϵ_q(q))
    Γ(i, j, p, q) = 3/2*Γ_i[1]*ϵ_q(p)*imag.(GreenTensor(positions[i] - positions[j])) * adjoint(ϵ_q(q))
    q_list = [-1, 0, 1]
    H           = sum(J(i, j, p, q)  * (dagger(Σ_iq(i, p, F_i[1], F_i[2]))*Σ_iq(j, q, F_i[1], F_i[2]))  for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list)
    Jump        = [Σ_iq(i, p, F_i[1], F_i[2])           for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    Jump_dagger = [dagger(Σ_iq(j, q, F_i[1], F_i[2]))  for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    rates       = [Γ(i, j, p, q)                  for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    indx       = [(i, j, p, q)                   for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    
    # Checking rate vector
    if check_rates
        println("Checking rates: ")
        for ii in eachindex(rates)
            println((indx[ii], real.(round.(rates[ii], digits=3))))
        end
        println("Checking Omega: ")
        for q in q_list
            for p in q_list
                for j in [1, 2]
                    for i in [1, 2]
                        println((i, j, p, q), real.(round.(J(i, j, p, q), digits=3)))
                    end
                end
            end
        end
    end

    p = Progress(length(tspan), "Solving master equation...")
    
    function forward_progress(t, rho) 
        # function for checking the progress
        next!(p)
        return copy(rho)
    end
    t_out, ρ_t = timeevolution.master(tspan, ρ_0, H, Jump, 
        Jdagger=Jump_dagger, 
        rates=rates, fout=forward_progress
    )
    result = copy(parameters)
    @pack! result = t_out, ρ_t, Jump, rates
    return result
end

"""
Plot time evolution data
"""
function plot_dynamics(result::Dict; kwargs...)
    @unpack ρ_t, t_out, F_i, positions, Γ_i= result
    # Plot population dynamics
    fig1 = plot(title="spin-e: $(latexify(F_i[1])), spin-g: $(latexify(F_i[2]))", xlab="tΓ",)

    rhoee_1 = [real.(
        expect( (identityoperator(SpinBasis(F_i[1]))    ⊕ 0*identityoperator(SpinBasis(F_i[2]))) ⊗ 
                (identityoperator(SpinBasis(F_i[1]))    ⊕ identityoperator(SpinBasis(F_i[2])))     , ρ_t[ii] ) 
            ) for ii in eachindex(ρ_t)]
    rhogg_1 =  [real.(
        expect( (0*identityoperator(SpinBasis(F_i[1]))    ⊕ identityoperator(SpinBasis(F_i[2]))) ⊗ 
                (identityoperator(SpinBasis(F_i[1]))    ⊕ identityoperator(SpinBasis(F_i[2])))     , ρ_t[ii] ) 
            ) for ii in eachindex(ρ_t)]
    rhoee_2 = [real.(
        expect( (identityoperator(SpinBasis(F_i[1]))    ⊕ identityoperator(SpinBasis(F_i[2]))) ⊗ 
                (identityoperator(SpinBasis(F_i[1]))    ⊕ 0*identityoperator(SpinBasis(F_i[2])))     , ρ_t[ii] ) 
            ) for ii in eachindex(ρ_t)]
    rhogg_2 =  [real.(
        expect( (identityoperator(SpinBasis(F_i[1]))    ⊕ identityoperator(SpinBasis(F_i[2]))) ⊗ 
                (0*identityoperator(SpinBasis(F_i[1]))    ⊕ identityoperator(SpinBasis(F_i[2])))     , ρ_t[ii] ) 
            ) for ii in eachindex(ρ_t)]
    plot!(fig1, t_out*Γ_i[1], rhoee_1, label="ptl1, e", ls=:solid, lw=2)
    plot!(fig1, t_out*Γ_i[1], rhogg_1, label="ptl1, g", ls=:solid, lw=2)
    plot!(fig1, t_out*Γ_i[1], rhoee_2, label="ptl2, e", ls=:solid, lw=2)
    plot!(fig1, t_out*Γ_i[1], rhogg_2, label="ptl2, g", ls=:solid, lw=2)
    plot!(fig1, t_out*Γ_i[1], (rhoee_2 + rhoee_1)/2, label="tot, e", ls=:dash, lw=2)
    plot!(fig1, t_out*Γ_i[1], (rhogg_2 + rhogg_1)/2, label="tot, g", ls=:dash, lw=2)
    plot!(fig1, t_out*Γ_i[1], real.(tr.(ρ_t)), label="Tr(ρ)", ls=:dash, lc=:black)
    plot!(leg=:outerright)
    # Position of the atoms
    fig2 = plot(xlabel="x", ylabel="y", zlabel="z", title="Atom positions", aspect_ratio=:equal)
    for ii in eachindex(positions)
        plot!(fig2, [positions[ii][1]], [positions[ii][2]], [positions[ii][3]], st=:scatter, label="Atom $(ii)")
    end
    fig = plot(fig1, fig2, layout=(2, 1), size=(600, 800))
    return fig
end