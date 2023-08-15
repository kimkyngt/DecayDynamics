### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ b4668ef1-9d19-4a8e-a6a1-262d649a55b7
using DrWatson

# ╔═╡ d35c4f2c-3b3d-11ee-2a72-f36026b81822
begin
	# Find and go to the project directory
	@quickactivate "DecayDynamics"
	projectname()
	cd(projectdir())
	pwd()
end

# ╔═╡ 71b3a8ce-fae4-4028-8262-619502add5c0
begin
	# Loading codes
	include(srcdir("Base.jl"))
	include(srcdir("single_atom_four_level.jl"))

	# Add packages
	Plots.default(fontfamily="Helvetica")
	md"Loading codes..."
end

# ╔═╡ b9a46890-5f35-46f7-a6f4-cb9b361d98e3
md"# B field strength check
Check the field strength without any decays
"

# ╔═╡ 25fcbbcb-ab64-4576-bba9-21b3e95476e6
function check_Lamor_precession()
	# Parameters for the simulation
	F_i = [1//2, 1//2, 1//2, 1//2] # two intermediate states, 3D1, 3P1, 3P0, 1S0
	g_i = [1, 1, 1, 1]
	Γ_kl = sparse(
	        [1, 1, 2],
	        [2, 3, 4],
	        [0*0.385, 0*0.597, 0.1], # Branching ratio, 3P1 assumed to be just 10 times slower. 
	        4, 4
	)
	frac = 1
	exc_frac = [frac, 1-frac] # Excitation fraction of level 1 and 3
	Bfield = [0, 0, 0.1]
	tspan = range(0, 20, length=300)
	ρ_0 = normalize(directsum(
	    sqrt(exc_frac[1])*(Fm_state(F_i[1], 1//2) + Fm_state(F_i[1], -1//2)) ,
	    Ket(SpinBasis(F_i[2])),
	    sqrt(exc_frac[2])*spinup(SpinBasis(F_i[3])),
	    Ket(SpinBasis(F_i[4]))
	))

	# Pack the parameters
	parameters = @dict F_i g_i Bfield Γ_kl ρ_0 tspan exc_frac

	# Solve the master equation
	result = evolve_master(parameters)
	figs_population = plot_dynamics(result)
	return plot(figs_population..., size=(1000, 1000), layout=(3, 3), margins=3Plots.mm)
end

# ╔═╡ 60884e7f-934a-4a35-a305-990aa328178e
check_Lamor_precession()

# ╔═╡ Cell order:
# ╠═b4668ef1-9d19-4a8e-a6a1-262d649a55b7
# ╟─d35c4f2c-3b3d-11ee-2a72-f36026b81822
# ╟─71b3a8ce-fae4-4028-8262-619502add5c0
# ╠═b9a46890-5f35-46f7-a6f4-cb9b361d98e3
# ╠═25fcbbcb-ab64-4576-bba9-21b3e95476e6
# ╠═60884e7f-934a-4a35-a305-990aa328178e
