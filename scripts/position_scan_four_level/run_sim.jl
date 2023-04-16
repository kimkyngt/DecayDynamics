using DrWatson
@quickactivate "DecayDynamics"
include("../../src/Base.jl")
include("../../src/many_atoms_four_level.jl")

# Folders to save data
foldername = "two_atoms_four_level"
subfoldername = "F_i=[2.5,1.5,1.5,0.5]m=2.5position_scan"

filecount = 1
m_exc = 5//2
@showprogress for _z in range(0.05, 3, length=101)
        global filecount, m_exc
        # Simulation parameters
        F_i = [5//2, 3//2, 3//2, 1//2] # two intermediate states, 3D1, 3P1, 3P0, 1S0
        Γ_kl = sparse(
                [1, 1, 2],
                [2, 3, 4],
                [1*0.385, 1*0.597, 0.1], # Branching ratio
                4, 4
        )
        knorm = 2π*sparse(
                [1, 1, 2],
                [2, 3, 4],
                [1, 2.6/2.7, 2.6/0.689], # wavevector ratio
                4, 4
        ) # length of the k-vector, reference to 2.6. index dictate which transition.
        frac = 1
        exc_frac = [frac, 1-frac] # Excitation fraction of level 1 and 3
        positions = [[0, 0, _z/2], [0, 0, -_z/2]] # 1 = 2.6 um. 
        displacement = positions[1]-positions[2] # for saving
        tspan = range(0, 37, length=100)
        ψsingle_exc = ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
        ρ_0 = normalize(ψsingle_exc ⊗ ψsingle_exc)
        description="Fully_excited"
        parameters = @dict F_i Γ_kl positions displacement ρ_0 tspan m_exc knorm exc_frac description

        # Run simulation
        @time result = evolve_master(parameters, check_rates=false)

        # plot and save
        fig = plot_dynamics(result)
        filename = savename(
                "[$(filecount)]"*description,
                parameters,
                connector="|", 
                allowedtypes=[Real, Int, Vector{<:Number}, ],
                sort=true,
                val_to_string=rational2str,
                accesses=[:description, :F_i, :m_exc, :displacement, :exc_frac],
        )

        savefig(fig, plotsdir(foldername, subfoldername, filename*".pdf"))
        wsave(datadir(foldername, subfoldername, filename*".jld2"), result)
        display(fig)
        filecount += 1
end

for _x in range(0.05, 3, length=101)
        global filecount
        # Simulation parameters
        F_i = [5//2, 3//2, 3//2, 1//2] # two intermediate states, 3D1, 3P1, 3P0, 1S0
        Γ_kl = sparse(
                [1, 1, 2],
                [2, 3, 4],
                [1*0.385, 1*0.597, 0.1], # Branching ratio
                4, 4
        )
        knorm = 2π*sparse(
                [1, 1, 2],
                [2, 3, 4],
                [1, 2.6/2.7, 2.6/0.689], # wavevector ratio
                4, 4
        ) # length of the k-vector, reference to 2.6. index dictate which transition.
        frac = 1
        exc_frac = [frac, 1-frac] # Excitation fraction of level 1 and 3
        positions = [[_x/2, 0, 0], [-_x/2, 0, 0]] # 1 = 2.6 um. 
        displacement = positions[1]-positions[2] # for saving
        tspan = range(0, 37, length=100)
        ψsingle_exc = ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
        ρ_0 = normalize(ψsingle_exc ⊗ ψsingle_exc)
        description="Fully_excited"
        parameters = @dict F_i Γ_kl positions displacement ρ_0 tspan m_exc knorm exc_frac description

        # Run simulation
        @time result = evolve_master(parameters, check_rates=false)

        # plot and save
        fig = plot_dynamics(result)
        filename = savename(
                "[$(filecount)]"*description,
                parameters,
                connector="|", 
                allowedtypes=[Real, Int, Vector{<:Number}, ],
                sort=true,
                val_to_string=rational2str,
                accesses=[:description, :F_i, :m_exc, :displacement, :exc_frac],
        )

        savefig(fig, plotsdir(foldername, subfoldername, filename*".pdf"))
        wsave(datadir(foldername, subfoldername, filename*".jld2"), result)
        display(fig)
        filecount += 1
end