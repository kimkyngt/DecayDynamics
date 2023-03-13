using Plots, QuantumOptics, CollectiveSpins, LaTeXStrings, WignerSymbols, ProgressMeter

""" 
Get |F, m>
"""
function Fm_state(F, m)
    b = SpinBasis(F)
    return normalize(sigmam(b)^Int(F - m) *spinup(b))
end

""" 
Atomic coherence operator.
"""
function σ(F_1, m_1, F_2, m_2)
    ψe = Fm_state(F_1, m_1)⊕ Ket(SpinBasis(F_2))
    ψg = Ket(SpinBasis(F_1)) ⊕ Fm_state(F_2, m_2)
    return sparse(ψe ⊗ dagger(ψg))
end

"""
Atomic lowering operator with hyperfine structure. Spherical basis.
"""
function Σ_q(q::Int, F_e, F_g)
    # Prevent the case where q-m is larger than F_e
    return sparse(dagger(sum(clebschgordan(F_g, -m, 1, q, F_e)*σ(F_e, m - q, F_g, m) for m = -F_g:1:F_g if abs(q-m) <= F_e )))
end

"""
Polarization vector of the light.
"""
function ϵ_q(q::Int)
    if q == 0
        return [0, 0, 1]
    elseif abs(q) == 1
        return -q/sqrt(2) * [1, q*im, 0]
    else 
        throw(ArgumentError("Argument q must be one of [-1, 0, 1]"))
    end
end

"""
Electric field operator.
"""
function E⁺(r::Vector, F_e, F_g)
    if norm(r) == 0
        G = sparse(zeros(3, 3) + 1im*ones(3, 3))
    else
        G = sparse(GreenTensor(r))
    end
    Gep(q) = G*conj(ϵ_q(q))
    return [sum(Gep(q)[ii] * Σ_q(q, F_e, F_g) for q in [-1, 0, 1]) for ii in 1:3]
end

"""
Intensity at r
"""
function get_intensity(r::Vector, ρ::Operator)
    F_e = ρ.basis_l.bases[1].spinnumber
    F_g = ρ.basis_l.bases[2].spinnumber
    real(expect(E⁺(r, F_e, F_g)' * E⁺(r, F_e, F_g), ρ))
end

"""
Draw radiation pattern
"""
function draw_radiation_pattern(ρ::Operator, nang::Int=20,; kwargs...)
    θ = range(0, π, length = nang)
    ϕ = range(0, 2π, length  = nang)
    x, y, z = polar_2_xyz(θ, ϕ, 1.5*ones(length(θ)))

    intensity_r = zero(x)
    for I in CartesianIndices(x)
        i, j = I[1], I[2]
        intensity_r[i, j] = get_intensity([x[i, j], y[i, j], z[i, j]], ρ)
    end
    intensity_r = intensity_r/maximum(intensity_r) # normalize
    data_x, data_y, data_z = polar_2_xyz(θ, ϕ, intensity_r)

    plt = surface(
        data_x, data_y, data_z,
        color = palette(:tab10)[1],
        xlabel="x", ylabel="y", zlabel="z",
        colorbar=false,
        xlim=(-1.1, 1.1), 
        ylim=(-1.1, 1.1), 
        zlim=(-1.1, 1.1), 
        xticks=[-1, 0, 1],
        yticks=[-1, 0, 1],
        zticks=[-1, 0, 1],
        aspect_ratio=1,
        dpi=200,
        frame_style=:box,
        ;kwargs...
    )
    return plt
end


"""
Get animation using Plot.jl's @animate macro. 
"""
function get_animation(experiment::Dict;kwargs...)
    @unpack t_out, ρ_t = experiment
    p = Progress(length(t_out), "Generating gif...")
    anim = @animate for ii in eachindex(t_out)
        # Update the plot with the i-th frame
        tlabel = t_out[ii]
        draw_radiation_pattern(ρ_t[ii], title="Time = $tlabel/Γ")
        next!(p)
    end
    return anim
end


"""
Convert polar coordinate to Cartesian for plotting.
"""
function polar_2_xyz(θ, ϕ, r=1)
    x = [sin(θi)*cos(ϕj) for θi ∈ θ,  ϕj ∈ ϕ].*r
    y = [sin(θi)*sin(ϕj) for θi ∈ θ,  ϕj ∈ ϕ].*r
    z = [cos(θi) 	     for θi ∈ θ,  ϕj ∈ ϕ].*r
    return x, y, z
end


"""
Solve time dynamics with master equation.
"""
function evolve_master(parameters::Dict)
    @unpack F_e, F_g, tspan, field_z, field_x, field_y, g_e, g_g, ρ_0 = parameters
    b_e = SpinBasis(F_e)
    b_g = SpinBasis(F_g)
    q_list = [-1, 0, 1]
    H = sum(interaction.Omega([0, 0, 0], [0, 0, 0], ϵ_q(q), ϵ_q(q)) * dagger(Σ_q(q, F_e, F_g)) * Σ_q(q, F_e, F_g) for q in q_list)
    H += 2*π*field_x*(g_e*sigmax(b_e)  ⊕ g_g*sigmax(b_g))/2
    H += 2*π*field_y*(g_e*sigmay(b_e)  ⊕ g_g*sigmay(b_g))/2
    H += 2*π*field_z*(g_e*sigmaz(b_e)  ⊕ g_g*sigmaz(b_g))/2
    H = sparse(H)
    J = [sqrt(interaction.Gamma([0, 0, 0], [0, 0, 0], ϵ_q(q), ϵ_q(q)))*Σ_q(q, F_e, F_g) for q in q_list]
    t_out, ρ_t = timeevolution.master_h(tspan, ρ_0, H, J)
    result = copy(parameters)
    @pack! result = t_out, ρ_t
    return result
end

"""
Plot time evolution data
"""
function plot_dynamics(parameters::Dict; kwargs...)
    @unpack ρ_t, t_out, F_e, F_g = parameters
    # Plot population dynamics
    palette_e = cgrad(:jet, Int(2*F_e+1), categorical=true)
    palette_g = cgrad(:jet, Int(2*F_g+1), categorical=true)
    b_e = SpinBasis(F_e)
    b_g = SpinBasis(F_g)
    fig1 = plot(ylab="Excited, $F_e")
    plot!(fig1, t_out, real.(expect(one(b_e) ⊕ 0*one(b_g), ρ_t)), label="Tr(ρᵉᵉ)", legend=:right, ls=:dash, lc=:black)
    for ii in range(1, Int(F_e*2+1))
        spin = F_e - (ii-1)
        state_to_plot = dagger(normalize((sigmam(b_e)^(ii-1) * spinup(b_e)) ⊕ Ket(b_g) ))
        plot!(fig1, t_out, real.(expect(projector(state_to_plot), ρ_t)), label="m = $spin", color=palette_e[ii])
    end
    fig2 = plot(ylab="Ground, $F_g", xlab="Time", leg=:right)
    plot!(fig2, t_out, real.(expect(0*one(b_e) ⊕ one(b_g), ρ_t)), label="Tr(ρᵍᵍ)", legend=:right, ls=:dash, lc=:black)
    for ii in range(1, Int(F_g*2+1))
        spin = F_g - (ii-1)
        state_to_plot = dagger(normalize( Ket(b_e) ⊕ (sigmam(b_g)^(ii-1) * spinup(b_g)) ))
        plot!(fig2, t_out, real.(expect(projector(state_to_plot), ρ_t)), label="m = $spin", color=palette_g[ii])
    end
    # excited, coherene
    fig3 = plot(ylab="Coherence, e", xlab="Time", leg=:right)
    for ii in range(1, Int(F_e*2))
        spin_l = F_e - (ii-1)
        spin_h = F_e - ii
        state_to_plot = (Fm_state(F_e, spin_h) ⊗ dagger(Fm_state(F_e, spin_l))) ⊕ Operator(b_g, b_g, zeros(length(b_g), length(b_g))) 
        plot!(fig3, t_out, real.(expect(state_to_plot, ρ_t)), label="ρᵉᵉ($spin_h, $spin_l)", color=palette_e[ii])
    end

    fig4 = plot(ylab="Coherence, g", xlab="Time", leg=:right)
    for ii in range(1, Int(F_g*2))
        spin_l = F_g - (ii-1)
        spin_h = F_g - ii
        state_to_plot = Operator(b_e, b_e, zeros(length(b_e), length(b_e))) ⊕ (Fm_state(F_g, spin_h) ⊗ dagger(Fm_state(F_g, spin_l)))
        plot!(fig4, t_out, real.(expect(state_to_plot, ρ_t)), label="ρᵍᵍ($spin_h, $spin_l)", color=palette_g[ii])
    end

    return plot(fig1, fig2, fig3, fig4, layout=(2, 2), size=(1200, 800), margins=5Plots.mm;kwargs...)
end