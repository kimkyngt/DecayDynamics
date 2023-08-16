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
    H  = directsum([Bfield[1]*g_i[ii]*1/2*sigmax(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
    H += directsum([Bfield[2]*g_i[ii]*1/2*sigmay(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
    H += directsum([Bfield[3]*g_i[ii]*1/2*sigmaz(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
    H = sparse(H)

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
    level_label = [L"^3D_1", L"^3P_1", L"^3P_0", L"^1S_0"]
    
    figs = []
    # Plot population
    for indx in eachindex(F_i)
        fig = plot(ylab="Population, $(level_label[indx]), $(F_i[indx])", legend=:best, xlab="tΓ")
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
        fig = plot(ylab="Coherence, $(level_label[indx]), $(F_i[indx])", legend=:best, xlab="tΓ")
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


"""
Electric field generate by F_i[2] to F_i[4] transition.
"""
function E⁺(r::Vector, F_i::Vector{<:Rational})
    if norm(r) == 0
        G = sparse(zeros(3, 3) + 1im*ones(3, 3))
    else
        G = sparse(GreenTensor(r))
    end
    Gep(q) = G* ϵ_q(q)'
    return [sum(Gep(q)[ii]*Σ_q(q, F_i, 2, 4) for q in [-1, 0, 1]) for ii in 1:3]
end

"""
<ψ|E⁽⁻⁾⋅E⁽⁺⁾|ψ>
"""
function get_intensity(r::Vector, ρ::Operator, F_i::Vector{<:Rational})
    real(expect(E⁺(r, F_i)' * E⁺(r, F_i), ρ))
end

function draw_radiation_pattern_3D(ρ::Operator, F_i::Vector{<:Rational}, nang::Int=20,; radius::Number=1e3, kwargs...)
    θ = range(0, π, length = nang)
    ϕ = range(0, 2π, length  = nang)
    x, y, z = polar_2_xyz(θ, ϕ, radius*ones(length(θ)))

    intensity_r = zero(x)
    for I in CartesianIndices(x)
        i, j = I[1], I[2]
        intensity_r[i, j] = get_intensity([x[i, j], y[i, j], z[i, j]], ρ, F_i)
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

function draw_radiation_pattern_2D(ρ::Operator, F_i::Vector{<:Rational}, nang::Int=20,; radius::Number=1, kwargs...)
    θ = range(0, π, length = nang)
    ϕ = range(0, 2π, length  = nang)
    x, y, z = polar_2_xyz(θ, ϕ, radius*ones(length(θ)))

    intensity_r = zero(x)
    for I in CartesianIndices(x)
        i, j = I[1], I[2]
        intensity_r[i, j] = get_intensity([x[i, j], y[i, j], z[i, j]], ρ, F_i)
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
Average the detector signal over detector NA. 
"""
function plot_detector_signal(result::Dict, r::Number, θ_center::Number, ϕ_center::Number; NA::Number=0.2, Nsample::Int=10)
    @unpack F_i, ρ_t, t_out = result
    # Sample detector position
    x, y, z = sample_omega_on_sphere(Nsample, NA, θ_center, ϕ_center)
    # Plot detector signal
    It = []
    for jj in eachindex(ρ_t)
        push!(It, sum(get_intensity(r*[x[ii], y[ii], z[ii]], ρ_t[jj], F_i) for ii in eachindex(x))/Nsample)
    end
    fig_sig = plot(t_out, It, xlab="tΓ", ylab="Detector signal", legend=false, color=1)
    # Plot detector gemetry
    fig_det = draw_unit_sphere()
    ang_to_print = Int(round(θ_center/π*180))
    plot!(fig_det, x, y, z, st=:scatter, color=1, xlab="y", ylab="x", zlab="z", ticks=false, title="Detector: NA=$NA, θ=$(ang_to_print) deg")
    return [fig_sig, fig_det, It, t_out]
end


"""
Getting the detector signals
"""
function get_detector_signal(result::Dict, r::Number, θ_center::Number, ϕ_center::Number; NA::Number=0.2, Nsample::Int=10)
    @unpack F_i, ρ_t, t_out = result
    # Sample detector position
    x, y, z = sample_omega_on_sphere(Nsample, NA, θ_center, ϕ_center)
    # Plot detector signal
    It = []
    for jj in eachindex(ρ_t)
        push!(It, sum(get_intensity(r*[x[ii], y[ii], z[ii]], ρ_t[jj], F_i) for ii in eachindex(x))/Nsample)
    end
    return t_out, It
end

"""
Double exponential model
"""
function double_exp(t)
    return (exp(-t/10) - exp(-t))
end
