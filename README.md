# DecayDynamics

## Single particle, multi level
- Notation following [2].

$$\begin{gathered}
\mathscr{H}=\hbar \sum_{q = -1}^1 J_{q} \hat{\Sigma}_{q}^{\dagger} \hat{\Sigma}_{ q^{\prime}}, \\
\mathscr{L}[\rho]=\sum_{q = -1}^1 \frac{\Gamma_{q}}{2}\left(2 \hat{\Sigma}_{ q^{\prime}} \rho \hat{\Sigma}_{q}^{\dagger}-\hat{\Sigma}_{q}^{\dagger} \hat{\Sigma}_{j q^{\prime}} \rho -\rho \hat{\Sigma}_{q}^{\dagger} \hat{\Sigma} \right) \\
J_{q}=0 \\
\Gamma_{q q}=\frac{2 \mu_0 \omega_0^2}{\hbar}|\wp|^2 \hat{\mathbf{e}}_q \cdot \operatorname{Im} \mathbf{G}\left(\mathbf{0}, \mathbf{0}, \omega_0\right) \cdot \hat{\mathbf{e}}_{q}^* \\
\mathbf{G}\left(\mathbf{r}, \omega_0\right)=\frac{e^{\mathrm{i} k_0 r}}{4 \pi k_0^2 r^3}\left[\left(k_0^2 r^2+\mathrm{i} k_0 r-1\right) \mathbb{1} + \left(-k_0^2 r^2-3 \mathrm{i} k_0 r+3\right) \frac{\mathbf{r} \otimes \mathbf{r}}{r^2}\right] \\
\hat{\Sigma}_{q}^{\dagger}=\sum_{m_g=-F_g}^{F_g} C_{m_g, q} \hat{\sigma}_{F_e m_g-q, F_g m_g} \\
\hat{\sigma}_{F_e m_g-q, F_g m_g}=\left|F_e m_g-q \right> \left< F_g m_g\right| \\
C_{m_g, q}=(-1)^{F_g-m_g}\left(\begin{array}{llc}
F_g & 1 & F_e \\
-m_g & q & m_g-q
\end{array}\right)
\end{gathered}
$$

$C_{m_g, q}$ is the Clebsch–Gordan coefficient, $\hat{\mathbf{e}}_q$ is the polarization vector in spherical coordinate. 

To calculate the radiation power at $\mathbf{r}$, evaluate 
$$
I(\mathbf{r}) = \left< \psi \right| \hat{\mathbf{E}}^{-} (\mathbf{r}) \cdot \hat{\mathbf{E}}^{+} (\mathbf{r}) \left| \psi \right>
$$
, where
$$
\hat{\mathbf{E}}^{+}(\mathbf{r})=\mu_0 \omega_0^2 \sum_{q=-1}^1 \mathbf{G}\left(\mathbf{r}, 0, \omega_0\right) \cdot \hat{\mathbf{e}}_q^* \wp \hat{\Sigma}_{j q} 
$$
.

## Cascaded multilevel
For 3D1-3P1-1S0 decay.

<!-- 
## Many particle, multi level
- Notation following [2].

$$\begin{gathered}
\mathscr{H}=\hbar \sum_{i, j=1}^N \sum_{q, q^{\prime}=-1}^1 J_{i j q q^{\prime}} \hat{\Sigma}_{i q}^{\dagger} \hat{\Sigma}_{j q^{\prime}}, \\
\mathscr{L}[\rho]=\sum_{i, j=1}^N \sum_{q, q^{\prime}=-1}^1 \frac{\Gamma_{i j q q^{\prime}}}{2}\left(2 \hat{\Sigma}_{j q^{\prime}} \rho \hat{\Sigma}_{i q}^{\dagger}-\hat{\Sigma}_{i q}^{\dagger} \hat{\Sigma}_{j q^{\prime}} \rho\right. \\
\left.-\rho \hat{\Sigma}_{i q}^{\dagger} \hat{\Sigma}_{j q^{\prime}}\right)
\end{gathered}$$
$$
\begin{aligned}
& J_{i j q q^{\prime}}=-\frac{\mu_0 \omega_0^2}{\hbar}|\wp|^2 \hat{\mathbf{e}}_q \cdot \operatorname{Re} \mathbf{G}\left(\mathbf{r}_i, \mathbf{r}_j, \omega_0\right) \cdot \hat{\mathbf{e}}_{q^{\prime}}^*, \\
& \Gamma_{i j q q^{\prime}}=\frac{2 \mu_0 \omega_0^2}{\hbar}|\wp|^2 \hat{\mathbf{e}}_q \cdot \operatorname{Im} \mathbf{G}\left(\mathbf{r}_i, \mathbf{r}_j, \omega_0\right) \cdot \hat{\mathbf{e}}_{q^{\prime}}^*,
\end{aligned}
$$
$$
\begin{aligned}
& \mathbf{G}\left(\mathbf{r}, \omega_0\right)=\frac{e^{\mathrm{i} k_0 r}}{4 \pi k_0^2 r^3}\left[\left(k_0^2 r^2+\mathrm{i} k_0 r-1\right) \mathbb{1}+\left(-k_0^2 r^2-3 \mathrm{i} k_0 r+3\right) \frac{\mathbf{r} \otimes \mathbf{r}}{r^2}\right], \\
&
\end{aligned}
$$

To calculate the electric field, 
$$
\hat{\mathbf{E}}^{+}(\mathbf{r})=\mu_0 \omega_0^2 \sum_{j=1}^N \sum_{q=-1}^1 \mathbf{G}\left(\mathbf{r}, \mathbf{r}_j, \omega_0\right) \cdot \hat{\mathbf{e}}_q^* \wp \hat{\Sigma}_{j q}
$$ -->

## Reference
1.  D. A. Steck, Quantum and Atom Optics (2022).
2.  A. Asenjo-Garcia, H. J. Kimble, D. E. Chang, Proc. Natl. Acad. Sci. U.S.A. 116, 25503–25511 (2019).
3.  B. Zhu, J. Cooper, J. Ye, A. M. Rey, Phys. Rev. A. 94, 023612 (2016).


<!-- This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> DecayDynamics

It is authored by Kyungtae Kim.

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
