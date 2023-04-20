# DecayDynamics

Contents moved to `docs/hyperfine_state_contribution.tex`.

## Reference

1. D. A. Steck, Quantum and Atom Optics (2022).

2. A. Asenjo-Garcia, H. J. Kimble, D. E. Chang, Proc. Natl. Acad. Sci. U.S.A. 116, 25503–25511 (2019).

3. B. Zhu, J. Cooper, J. Ye, A. M. Rey, Phys. Rev. A. 94, 023612 (2016).

4. [`CollectiveSpins.jl`](https://qojulia.github.io/CollectiveSpins.jl/dev/descriptions/) theoretical description page.

## Activating this project

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:

```julia
julia> using Pkg
julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
julia> Pkg.activate("path/to/this/project")
julia> Pkg.instantiate()
```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:

```julia
using DrWatson
@quickactivate "DecayDynamics"
```

which auto-activate the project and enable local path handling from DrWatson.
<!-- This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> DecayDynamics

It is authored by Kyungtae Kim.
<!-- 
To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "DecayDynamics"
```
which auto-activate the project and enable local path handling from DrWatson. -->
