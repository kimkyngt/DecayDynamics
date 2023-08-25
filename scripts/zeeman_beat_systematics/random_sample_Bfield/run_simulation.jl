# This script is for generating data for the systematics estimation based on the test results in notebooks/Zeeman_beat_test.ipynb
using DrWatson

# Loading codes
include(srcdir("single_atom_four_level.jl"))
include(srcdir("Base.jl"))


# Fixed parameters
F_i = [11//2, 9//2, 9//2, 9//2] # two intermediate states, 3D1, 3P1, 3P0, 1S0
g_i = [127, 85, 0, 0]/76 * 2π # 2π from the angular frequency - frequency conversion.
branching_ratio = [0.385, 0.597] / (0.385 + 0.597)
Γ_kl = sparse(
        [1, 1, 2],
        [2, 3, 4],
        [1*branching_ratio[1], 1*branching_ratio[2], 0.1], # Branching ratio, 3P1 assumed to be just 10 times slower. 
        4, 4
)

# Simulation codes

"""
        # Generate detector signal sample.
    - This time, branching ratio has to be summed to 1.
    - To do this fast, we use very small NA. and 1 sample.
    - The angles are detector location with respect to the quantization axis when we polarize the atom. Since we ramp down the field very quickly a after the polarization, this axis and the residual bias field direction would be rather different. The best way to sample the B field noise is fixing the initial population and the detector angle based on the moment we polarize the atoms, and then sample the B field in various directions. Our detector is in the x-y plane, and 4 degree deviated from the x-axis, and the quantization axis is along the x-axis when we polarize the atoms. θ = 0, ϕ = 4 degree.
    - so far the initial state is hard-coded.

    args: 
        Bfield::Vector: B field vector in the unit of Gauss.

"""
function generate_sample(Bfield::Vector, ρ_0, θ, ϕ)
    tspan = range(0, 77, length=300)
    # Pack the parameters
    parameters = @dict F_i g_i Bfield Γ_kl ρ_0 tspan

    # Solve the master equation
    result = evolve_master(parameters)
    t_out, It = get_detector_signal(result, 1e3, θ, ϕ, NA=0, Nsample=1)
    return t_out, It, result
end

@. model(x, p) = p[1]*(exp(-(x - p[4])/p[2]) - exp(-(x - p[4])/p[3])) + p[5] 

function sim_detector(Bfield, ρ_0, θ, ϕ; out_fig=false)
    t_out, It, _ = generate_sample(Bfield, ρ_0, θ, ϕ)
    It = It / maximum(It)
    # plot and fit the data to a double exponential
    p0 = [maximum(It), 10, 1, 0, 0]
    fit_raw = curve_fit(model, t_out, It, p0)
    fit_model = model(t_out, fit_raw.param)
    se = LsqFit.stderror(fit_raw)

    if out_fig
        fig = plot(t_out, It, lab="simulation")
        plot!(t_out, fit_model, lab="double exponential fit", ls=:dash)
        fig_res = plot(t_out, fit_raw.resid, lab="fit residuals")

        figtot = plot(fig, fig_res, layout=(2, 1), size=(800, 600), suptitle="τD = $(round(fit_raw.param[3], digits=4))($(round(se[3], digits=4))), τP = $(round(fit_raw.param[2], digits=3))($(round(se[2], digits=4))), \n (θ, |B|) = ($(round(θ, digits=3)), $(round(norm(Bfield), digits=3)))", titlefontsize=8)

        return figtot
    else
        return fit_raw, t_out, It
    end
end


"""
Main simulation running codes. 
"""
function run_simulation(d::Dict)
    @unpack rho0, Bfield, thetaDet, phiDet=d
    fit_raw, t_out, It = sim_detector(Bfield, rho0, thetaDet, phiDet)
    tosave = @strdict fit_raw t_out It
    return tosave
end

N_simulation = 100
thetaDet = 4*π/180
phiDet = 0
rho0 = normalize(directsum(
    Fm_state(F_i[1], 9//2) + Fm_state(F_i[1], 5//2) ,
    Ket(SpinBasis(F_i[2])),
    Ket(SpinBasis(F_i[3])),
    Ket(SpinBasis(F_i[4]))))

σB = 2e-3 # standard deviation of the gaussian

@showprogress for ii in 1:N_simulation
    # Sampling Bfield
    x, y, z = sample_sphere(1)
    Bfield = randn()*σB*[x, y, z]
    params = @strdict rho0 Bfield thetaDet phiDet σB
    result = run_simulation(params)
    data_to_save = @strdict result params
    save_file_name = savename(params, "jld2")
    safesave(datadir("zeeman_beat", "gaussian_sampling", save_file_name), data_to_save)
    println(save_file_name)
    # show progress
    # println("$(round(100*ii / N_simulation, digits=2))%")
end
