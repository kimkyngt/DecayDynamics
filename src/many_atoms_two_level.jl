using Plots, QuantumOptics, LaTeXStrings, WignerSymbols, ProgressMeter

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
Solve time dynamics w ith master equation.
"""
function evolve_master(parameters::Dict)
    @unpack F_i, g_i, Γ_i, Bfield, kr, ρ_0, tspan = parameters

    J(i, j, p, q) = real(Γ_i[1]*(- ϵ_q(p)*real.(GreenTensor(kr[i] - kr[j])) * adjoint(ϵ_q(q))))
    Γ(i, j, p, q) = real(Γ_i[1]*(2*ϵ_q(p)*imag.(GreenTensor(kr[i] - kr[j])) * adjoint(ϵ_q(q))))
    q_list = [-1, 0, 1]
    H           = sum(  J(i, j, p, q)  * (dagger(Σ_iq(i, p, F_i[1], F_i[2]))*Σ_iq(j, q, F_i[1], F_i[2])) for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list)
    Jump        = [Σ_iq(i, p, F_i[1], F_i[2])          for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    Jump_dagger = [dagger(Σ_iq(j, q, F_i[1], F_i[2]))  for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    rates       = [Γ(i, j, p, q) for i in [1, 2] for j in [1, 2] for p in q_list for q in q_list]
    p = Progress(length(tspan), "Solving master equation...")
    function fout(t, rho) 
        # function for checking the progress
        next!(p)
        return rho
    end
    t_out, ρ_t = timeevolution.master(tspan, ρ_0, H, Jump, Jdagger=Jump_dagger, rates=rates, fout=fout)
    # t_out, ρ_t = timeevolution.master_h(tspan, ρ_0, H, J)
    result = copy(parameters)
    @pack! result = t_out, ρ_t
    return result
end

"""
Plot time evolution data
"""
function plot_dynamics(result::Dict; kwargs...)
    @unpack ρ_t, t_out, F_i, kr= result
    # Plot population dynamics
    fig1 = plot(ylab="Excited, $(F_i[1])", xlab="t/Γ", )

    plot!(fig1, t_out, real.(tr.(ρ_t)), label="Tr(ρ)", legend=:right, ls=:solid, lc=:black)
    plot!(fig1, t_out, [real.(
        expect( (identityoperator(SpinBasis(F_i[1]))⊕projector(Ket(SpinBasis(F_i[2])))) ⊗ 
                (identityoperator(SpinBasis(F_i[1]))⊕projector(Ket(SpinBasis(F_i[2])))), ρ_t[ii] ) 
            ) for ii in eachindex(ρ_t)], label="Tr(ρᵉᵉ)", legend=:right, ls=:dash, lc=:red
        )
    plot!(fig1, t_out, [real.(
        expect( (projector(Ket(SpinBasis(F_i[1]))) ⊕ identityoperator(SpinBasis(F_i[2]))) ⊗ 
                (projector(Ket(SpinBasis(F_i[1]))) ⊕ identityoperator(SpinBasis(F_i[2]))), ρ_t[ii] ) 
            ) for ii in eachindex(ρ_t)], label="Tr(ρᵍᵍ)", legend=:right, ls=:dash, lc=:blue
        )

    fig = fig1
    return fig
end