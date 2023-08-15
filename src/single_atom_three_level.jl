using Plots, QuantumOptics, CollectiveSpins, LaTeXStrings, WignerSymbols, ProgressMeter, LinearAlgebra
import Rotations: RotX, RotZ

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
function Σ_q(q::Int, F_1, F_2)
    return sparse(dagger(sum(clebschgordan(F_2, m, 1, q, F_1)*σ(F_1, m + q, F_2, m) for m = -F_2:1:F_2 if abs(m + q) <= F_1 )))
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
Electric field generate by F_i[2] to F_i[3] transition. For only angular dependence calculation.
"""
function E⁺(r::Vector, F_i::Vector{<:Rational})
    if norm(r) == 0
        G = sparse(zeros(3, 3) + 1im*ones(3, 3))
    else
        G = sparse(GreenTensor(r))
    end
    Gep(q) = G*conj(ϵ_q(q))
    return [sum(Gep(q)[ii]*(0*one(SpinBasis(F_i[1])) ⊕ Σ_q(q, F_i[2], F_i[3])) for q in [-1, 0, 1]) for ii in 1:3]
end

"""
<ψ|E⁽⁻⁾⋅E⁽⁺⁾|ψ>
"""
function get_intensity(r::Vector, ρ::Operator, F_i::Vector{<:Rational})
    real(expect(E⁺(r, F_i)' * E⁺(r, F_i), ρ))
end

"""
Draw radiation pattern.
"""
function draw_radiation_pattern(ρ::Operator, F_i::Vector{<:Rational}, nang::Int=20,; kwargs...)
    θ = range(0, π, length = nang)
    ϕ = range(0, 2π, length  = nang)
    x, y, z = polar_2_xyz(θ, ϕ, 1.5*ones(length(θ)))

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
Plot signal over all 4pi solid angle.
"""
function get_detector_signal(ρ::Operator, F_i::Vector{<:Rational}, ; kwargs...)
    θ = range(0, π, length = nang)
    ϕ = range(0, 2π, length  = nang)
    x, y, z = polar_2_xyz(θ, ϕ, ones(length(θ)))

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
Get animation using Plot.jl's @animate macro. 
"""
function get_animation(result::Dict;kwargs...)
    @unpack t_out, ρ_t, F_i= result
    p = Progress(length(t_out), "Generating gif...")
    anim = @animate for ii in eachindex(t_out)
        # Update the plot with the i-th frame
        tlabel = t_out[ii]
        draw_radiation_pattern(ρ_t[ii], F_i, title="Time = $tlabel/Γ")
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
Solve time dynamics w ith master equation.
"""
function evolve_master(parameters::Dict)
    @unpack F_i, g_i, Γ_i, Bfield, ρ_0, tspan = parameters
    q_list = [-1, 0, 1]
    H  = π*directsum([Bfield[1]*g_i[ii]*sigmax(SpinBasis(F_i[ii])) for ii in [1, 2, 3]]...)
    H += π*directsum([Bfield[2]*g_i[ii]*sigmay(SpinBasis(F_i[ii])) for ii in [1, 2, 3]]...)
    H += π*directsum([Bfield[3]*g_i[ii]*sigmaz(SpinBasis(F_i[ii])) for ii in [1, 2, 3]]...)
    H = sparse(H)
    J = vcat([sqrt(Γ_i[1])*(Σ_q(q, F_i[1], F_i[2]) ⊕ one(SpinBasis(F_i[3]))) for q in q_list], [sqrt(Γ_i[2])*(one(SpinBasis(F_i[1])) ⊕ Σ_q(q, F_i[2], F_i[3])) for q in q_list])
    p = Progress(length(tspan), "Solving master equation...")
    function fout(t, rho) 
        # function for checking the progress
        next!(p)
        return rho
    end
    # t_out, ρ_t = timeevolution.master_h(tspan, ρ_0, H, J, fout=fout)
    t_out, ρ_t = timeevolution.master_h(tspan, ρ_0, H, J)
    result = copy(parameters)
    @pack! result = t_out, ρ_t
    return result
end

"""
Plot time evolution data
"""
function plot_dynamics(result::Dict; kwargs...)
    @unpack F_i, g_i, Γ_i, Bfield, ρ_t, t_out = result

    b_i = [SpinBasis(F_i[ii]) for ii in [1, 2, 3]]
    
    figs = []
    # Plot population
    for indx in [1, 2, 3]
        fig = plot(ylab="Population, Level-$indx, $(F_i[indx])", legend=:best, xlab="t/Γ")
        id_ops = [one(b_i[ii]) for ii in [1, 2, 3]]
        id_vec = circshift([1, 0, 0], indx-1)
        # Plot population of level-indx
        plot!(fig, t_out, real.(expect(directsum([one(b_i[ii]) for ii in [1, 2, 3]]...), ρ_t)), label="Tr(ρ)", ls=:solid, lc=:black)
        plot!(fig, t_out, real.(expect(directsum([id_vec[ii]*id_ops[ii] for ii in [1, 2, 3]]...), ρ_t)), label="Tr(ρ$indx)",  ls=:dash, lc=:black)

        # Plot each sublevel's population
        mycolor = cgrad(:Spectral, Int(2*F_i[indx]+1), categorical=true)
        blank_state = [Ket(b_i[ii]) for ii in [1, 2, 3]]
        for ii in range(1, Int(F_i[indx]*2+1))
            spin = F_i[indx] - (ii-1)
            _blank_state = copy(blank_state)
            _blank_state[indx] = Fm_state(F_i[indx], spin)
            state_to_plot = directsum(_blank_state...)
            plot!(fig, t_out, real.(expect(projector(state_to_plot), ρ_t)), label="m = $spin", color=mycolor[ii])
        end
        push!(figs, fig)
    end
    for indx in [1, 2, 3]
        fig = plot(ylab="Coherence, Level-$indx, $(F_i[indx])", legend=:best, xlab="t/Γ")
        # Plot each sublevel's coherence
        # mycolor = cgrad(:Spectral, Int(2*F_i[indx]+1), categorical=true)
        blank_state = [projector(Ket(b_i[ii])) for ii in [1, 2, 3]]
        for ii in range(1, Int(3))
            spin_h = F_i[indx] - (ii-1)
            spin_l = F_i[indx] - ii
            spin_ll = F_i[indx] - (ii+1)
            _blank_state = copy(blank_state)
            _blank_state[indx] = projector(Fm_state(F_i[indx], spin_h), dagger(Fm_state(F_i[indx], spin_l)))
            state_to_plot = directsum(_blank_state...)
            plot!(fig, t_out, real.(expect(state_to_plot, ρ_t)), label="($spin_h, $spin_l)", 
            # color=mycolor[ii], 
            ylim=(-0.5, 0.5),)
            # plot!(fig, t_out, imag.(expect(state_to_plot, ρ_t)), label="", color=mycolor[ii], ls=:dash,  )
            _blank_state = copy(blank_state)
            _blank_state[indx] = projector(Fm_state(F_i[indx], spin_h), dagger(Fm_state(F_i[indx], spin_ll)))
            state_to_plot = directsum(_blank_state...)
            plot!(fig, t_out, real.(expect(state_to_plot, ρ_t)), label="($spin_h, $spin_ll)", 
            # color=mycolor[ii], 
            ls=:dot, ylim=(-0.5, 0.5),)
            # plot!(fig, t_out, imag.(expect(state_to_plot, ρ_t)), label="", color=mycolor[ii],  ls=:dashdot)
        end
        push!(figs, fig)
    end

    # return plot(figs..., size=(1500, 800), layout=(2, 3), margins=30Plots.px;kwargs...)
    return figs
end


"""
Plot time-dependent radiation power at (θ, ϕ)
"""
function get_detector_signal(result::Dict, θ::Number, ϕ::Number)
    @unpack F_i, ρ_t, t_out = result
    r = [cos(θ)*cos(ϕ), cos(θ)*sin(ϕ), sin(θ)]
    I_t = [get_intensity(r, ρ, F_i) for ρ in ρ_t]
    return t_out, I_t
end

"""
Randomly sample unit vectors in NA centered at (theta, phi).
Returns unit vector components x, y, z
test: plot(sample_omega_on_sphere(10, 0.5, pi/2, pi), st:scatter)
"""
function sample_omega_on_sphere(Nsample::Int, NA::Number, theta::Number, phi::Number)
    ang = asin(NA)
    theta_sample = 2*ang*(rand(Nsample) .- 0.5)
    phi_sample = 2*pi*rand(Nsample)
    u = []
    for ii in eachindex(theta_sample)
        push!(u, RotZ(phi)*RotX(theta)*[
            sin(theta_sample[ii])*cos(phi_sample[ii]), 
            sin(theta_sample[ii])*sin(phi_sample[ii]), 
            cos(theta_sample[ii])
            ])
    end
    x, y, z = [], [], []
    for ii in eachindex(u)
        push!(x, u[ii][1]) 
        push!(y, u[ii][2])
        push!(z, u[ii][3])
    end
    return x, y, z
end

"""
Draw a unit sphere. Nsample determines the density
"""
function draw_unit_sphere(Nsample=20, )
    θ = range(0, π, length = Nsample)
    ϕ = range(0, 2π, length  = Nsample)
    fig = plot(
        polar_2_xyz(θ, ϕ, ones(length(θ))),
        st=:surface,
        color=:grey, 
        legend=false,
        )
    return fig
end


"""
Average the detector signal over detector NA. 
"""
function plot_detector_signal(result::Dict, θ_center::Number, ϕ_center::Number, NA::Number=0.2, Nsample::Int=10;kwargs...)
    @unpack F_i, ρ_t, t_out = result
    # Sample detector position
    x, y, z = sample_omega_on_sphere(Nsample, NA, θ_center, ϕ_center)
    # Plot detector signal
    It = []
    for jj in eachindex(ρ_t)
        push!(It, sum(get_intensity([x[ii], y[ii], z[ii]], ρ_t[jj], F_i) for ii in eachindex(x))/Nsample)
    end
    fig_sig = plot(t_out, It, xlab="t/Γ", ylab="Detector signal", legend=false, color=1)
    # Plot detector gemetry
    fig_det = draw_unit_sphere()
    ang_to_print = Int(round(θ_center/π*180))
    plot!(fig_det, x, y, z, st=:scatter, color=1, xlab="y", ylab="x", zlab="z", ticks=false, title="Detector: NA=$NA, θ=$(ang_to_print)ᵒ")
    return [fig_sig, fig_det]
end


"""
Just for generating filename
"""
function strfloat(x)
    a = string.(round.(float.(x), digits=3))
    if a isa String
            return a
    else
            return join(a, ",")
    end
end