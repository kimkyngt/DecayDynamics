using DrWatson
@quickactivate "DecayDynamics"
include("../../src/Base.jl")
include("../../src/many_atoms_four_level.jl")

foldername = "two_atoms_four_level"
subfoldername = "fraction_scan_0.21S0"

# Computing collective effects of two atoms with three level hyperfine structure
F_i = [2//1, 1//1, 1//1, 0//1] # two intermediate states, 3D1, 3P1, 3P0, 1S0

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
fracs = range(0, 1, length=52)
filecount = 1
tauDratio_list = []
tauPratio_list = []
chisq_list = []
for frac in fracs
        global filecount
        exc_frac = 0.8*[frac, 1-frac] # Excitation fraction of level 1 and 3
        # Bfield = [0, 0, 0]
        # g_i = [1, 1, 1, 1]
        positions = [[0, 0, 00.00], [0, 0, 0]] # 1 = 2.6 um. 
        displacement = positions[1]-positions[2] # for saving
        tspan = range(0, 20, length=300)
        m_exc = F_i[1]  # excited state

        # ρ_0 = ( ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...) ⊗ 
        #         ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...) )
        # description="m=$(AbstractFloat(m_exc))_FullyExcitedState_r12=$(round.(positions[1]-positions[2], digits=2))"

        ψsingle_exc = sqrt(exc_frac[1])*spinup(SpinBasis(F_i[1])) ⊕ Ket(SpinBasis(F_i[2])) ⊕ sqrt(exc_frac[2])*spinup(SpinBasis(F_i[3])) ⊕ sqrt(0.2)*spinup(SpinBasis(F_i[4]))

        # ψsingle_exc = ⊕([ii==1 ? Fm_state(F_i[ii], m_exc) : Ket(SpinBasis(F_i[ii])) for ii in eachindex(F_i)]...)
        ρ_0 = normalize(ψsingle_exc ⊗ ψsingle_exc)
        # description="Fully_excited"

        parameters = @dict F_i Γ_kl positions displacement ρ_0 tspan m_exc knorm exc_frac
        @time result = evolve_master(parameters, check_rates=false)
        fig, fitresult = plot_dynamics(result, return_dict=true)
        function strfloat(x)
                a = string.(round.(float.(x), digits=3))
                if a isa String
                        return a
                else
                        return join(a, ",")
                end
        end
        filename = savename(
                "[$(filecount)]",
                parameters,
                connector="|", 
                allowedtypes=[Real, Int, Vector{<:Number}, ],
                sort=true,
                val_to_string=strfloat,
                accesses=[:m_exc, :F_i, :displacement, :exc_frac],
        )

        savefig(fig, plotsdir(foldername, subfoldername, filename*".pdf"))
        # wsave(datadir(foldername, subfoldername, filename*".jld2"), tostringdict(result))
        push!(tauDratio_list, fitresult["tauDratio"])
        push!(tauPratio_list, fitresult["tauPratio"])
        push!(chisq_list, fitresult["chisq"])
        println(filename)
        display(fig)
        filecount += 1
end

# fig = plot(fracs[2:end], tauDratio_list[2:end], xlab="Excitation fraction", ylab="τfit/τ3D1", legend=false)
# savefig(fig, plotsdir(foldername, subfoldername, "tauDratio.pdf"))

fig1 = plot(fracs[2:end], tauDratio_list[2:end], xlab="Excitation fraction (3D1/3P0)", ylab="τfit/τ_i",  size=(400, 300),lab="i=3D1")
plot!(fracs[2:end], tauPratio_list[2:end], lab="i=3P1", legend=:right,)
savefig(fig1, plotsdir(foldername, subfoldername, "tauratio.pdf"))
display(fig1)
