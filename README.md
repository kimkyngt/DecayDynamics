# DecayDynamics

## Single particle, multi level

- Notation following [2].

$$\begin{gathered}
\mathscr{H}= \sum_{q = -1}^1 \sum_{i = g, e} \mu_B g_i \mathbf{F}^i \cdot \mathbf{B} \\
\mathscr{L}[\rho]=\sum_{q = -1}^1 \frac{\Gamma_{q}}{2}\left(2 \hat{\Sigma}_{ q^{\prime}} \rho \hat{\Sigma}_{q}^{\dagger}-\hat{\Sigma}_{q}^{\dagger} \hat{\Sigma}_{q} \rho -\rho \hat{\Sigma}_{q}^{\dagger} \hat{\Sigma}_{q} \right) \\
J_{q}=0 \\
\Gamma_{q} = 1 \\
\hat{\Sigma}_{q}^{\dagger}=\sum_{m_g=-F_g}^{F_g} C_{m_g, q}^{F_e, F_g} \hat{\sigma}_{F_e m_g+q, F_g m_g} \\
\hat{\sigma}_{F_e m_g+q, F_g m_g}=\left|F_e m_g+q \right> \left< F_g m_g\right| \\
C_{m_g, q}^{F_e, F_g} = \left< F_g, m_g;1, q|F_e, m_g+q\right>
\end{gathered}$$

$C_{m_g, q}$ is the Clebsch–Gordan coefficient, $\hat{\mathbf{e}}_q$ is the polarization vector in Cartesian coordinate. $\Gamma_q$ does not depend on a specific choise of the hyperfine state but depends on the fine structure state.

To calculate the radiation power at $\mathbf{r}$, evaluate 

$$\begin{gathered}
I(\mathbf{r}) = \left< \psi \right| \hat{\mathbf{E}}^{-} (\mathbf{r}) \cdot \hat{\mathbf{E}}^{+} (\mathbf{r}) \left| \psi \right>
\end{gathered}$$

, where the electric field operator is defined as

$$\begin{gathered}
\hat{\mathbf{E}}^{+}(\mathbf{r})= \mu_{0} \omega_{0}^2 \sum_{q = -1}^1 \mathbf{G} \left(\mathbf{r}, 0, \omega_{0} \right) \cdot \hat{\mathbf{e}}_{q}^{*} \wp \hat{\Sigma}_{j q} \\
\mathbf{G}\left(\mathbf{r}, \omega_0\right)=\frac{e^{\mathrm{i} k_0 r}}{4 \pi k_0^2 r^3}\left[\left(k_0^2 r^2+\mathrm{i} k_0 r-1\right) \mathbb{1} + \left(-k_0^2 r^2-3 \mathrm{i} k_0 r+3\right) \frac{\mathbf{r} \otimes \mathbf{r}}{r^2}\right]
\end{gathered}$$


## Cascaded multilevel
### Schematics

```mermaid
flowchart TD
    A(3D1, F = F_1, level-1) -->|2 μs| B(3P1, F = F_2, level-2) --> |20 μs| C(1S0, F = F_3, level-3)
```

First we define atomic lowering operator again.

$$
\begin{gathered}
{\hat{\Sigma}_{q}^{F_1, F_2}}^{\dagger}=\sum_{m_2=-F_2}^{F_2} C_{m_2, q}^{F_1, F_2} \hat{\sigma}_{F_1, m_2+q, F_2, m_2} \\
C_{m_2, q}^{F_1, F_2} = \left< F_2, m_2;1, q|F_1, m_2+q\right>
\end{gathered}
$$

Green function is depending on the wavelength of the transition.

$$ 
\mathbf{G}\left(\mathbf{r}, \omega_i\right)=\frac{e^{\mathrm{i} k_i r}}{4 \pi k_i^2 r^3}\left[\left(k_i^2 r^2+\mathrm{i} k_i r-1\right) \mathbb{1} + \left(-k_i^2 r^2-3 \mathrm{i} k_i r+3\right) \frac{\mathbf{r} \otimes \mathbf{r}}{r^2}\right]
$$

However, this will not matter as long as we are considering a single atom. It will not change the angular dependence of the radiation pattern as well. 

We set the decay rate of level-$i$ $\Gamma_q^i$ to be 10(1) for the level-1(2).

$$
\mathscr{L}[\rho]=\sum_{i = 1}^2 \sum_{q = -1}^1 \frac{\Gamma_{q}^{i}}{2}\left(2 \hat{\Sigma}_{ q^{\prime}}^{F_i, F_{i+1}} \rho {\hat{{\Sigma}_{q}}^{F_i, F_{i+1}}}^{\dagger}-{\hat{\Sigma}_{q}^{F_i, F_{i+1}}}^{\dagger} \hat{\Sigma}_{q}^{F_i, F_{i+1}} \rho -\rho {\hat{\Sigma}_{q}^{F_i, F_{i+1}}}^{\dagger} \hat{\Sigma}_{q}^{F_i, F_{i+1}} \right)
$$

For the Hamiltonian, we have rotating terms.
$$
\mathscr{H}=\hbar \sum_{i = 1}^3 \mu_B g_i \mathbf{F}^i \cdot \mathbf{B}
$$

<!-- ### Laser excitation

We expect geometric constraints of the superposition state of level-1.  -->

## Branching ratio
From the reference [1], equation (7.283), the reduced dipole moment of the transition between two different fine structure states are given by


$$\begin{aligned}
\left< J\|\mathbf{d}\| J^{\prime}\right> & \equiv \left< L S J\|\mathbf{d}\| L^{\prime} S J^{\prime} \right> \\
& = \left< L\|\mathbf{d}\| L^{\prime} \right> (-1)^{J^{\prime}+L+1+S} \sqrt{\left( 2 J^{\prime}+1 \right) (2 L+1)} \left\{ \begin{array}{ccc}
L & L^{\prime} & 1 \\
J^{\prime} & J & S \\
\end{array} \right\}
\end{aligned}$$

Here, (un)primed numbers are for the (excited)ground states. For ${}^{3}{D}_{1} \rightarrow {}^3{P}_{J^{\prime}}, \quad J^{\prime} \in (0, 1, 2)$ decay, $J = 1$, $L=2$, $L'=1$, $S=1$ are fixed. The ratio between squared dipole matrix elements of different ${}^3{P}_{J^{\prime}}$ states will be 

$$
\left| \left< {}^{3}{D}_{1} \| \mathbf{d} \| {}^3{P}_{J^{\prime}} \right> \right|^2 = (2J^{\prime}+1)(2\cdot 2+1)\left| \left\{\begin{array}{ccc}
2 & 1 & 1 \\
J^{\prime} & 1 & 1 \\
\end{array}\right\}  \right|^2 
$$

This gives the ratio $\frac{5}{9}:\frac{5}{12}:\frac{1}{36}$. The actual decay rate will depends on the frequency of the transition. For the decay rate $\Gamma$ from $e$ to $g$, 

$$
\Gamma = \frac{\omega_0^3}{3\pi \epsilon_0 \hbar c^3}| \left< g|\mathbf{d}|e \right>|^2.
$$

The decay rate will be proportional to the cube of the frequency or the inverse cube of the wavelength. The transition wavelengths of ${}^3{P}_{J^{\prime}}$ are $(2.6, 2.74, 3.07)~\mu m$. The actual decay ratio, considering $\frac{1}{2.6^3}:\frac{1}{2.74^3}:\frac{1}{3.07^3} \approx 0.69:0.59:0.42$, will be

$$
\frac{5}{9}\frac{1}{2.6^3}:\frac{5}{12}\frac{1}{2.74^3}:\frac{1}{36}\frac{1}{3.07^3}
$$

, which results in decay rates very close to [2].

We now consider the hyperfine structure. The matrix element has the same form as before. 

$$
\begin{aligned}
\left< F\|\mathbf{d}\| F^{\prime}\right> & \equiv\left< J I F\|\mathbf{d}\| J^{\prime} I F^{\prime}\right> \\
& =(-1)^{F^{\prime}+J+1+I} \sqrt{\left(2 F^{\prime}+1\right)(2 J+1)}\left\{\begin{array}{ccc}
J & J^{\prime} & 1 \\
F^{\prime} & F & I
\end{array}\right\} \\
& = \left< L\|\mathbf{d}\| L^{\prime}\right>(-1)^{J^{\prime}+L+1+S} \sqrt{\left(2 J^{\prime}+1\right)(2 L+1)}\left\{\begin{array}{ccc}
L & L^{\prime} & 1 \\
J^{\prime} & J & S \\
\end{array}\right\} \\
& \times (-1)^{F^{\prime}+J+1+I} \sqrt{\left(2 F^{\prime}+1\right)(2 J+1)}\left\{\begin{array}{ccc}
J & J^{\prime} & 1 \\
F^{\prime} & F & I
\end{array}\right\} .
\end{aligned}
$$

$$
\begin{aligned}
\left| \left< {}^{3}{D}_{1}, F\|\mathbf{d}\| {}^3{P}_{J^{\prime}}, F^{\prime} \right> \right|^2 =& (2J^{\prime} + 1)(2 \cdot 2 + 1)(2F^{\prime} + 1)(2 \cdot 1 + 1) \\ 
& \times \left|
   \left\{\begin{array}{ccc}
2 & 1 & 1 \\
J^{\prime} & 1 & 1 \\
\end{array}\right\}
\left\{\begin{array}{ccc}
1 & J^{\prime} & 1 \\
F^{\prime} & F & 9/2
\end{array}\right\}\right|^2 \left| \left< L\|\mathbf{d}\| L^{\prime}\right> \right|^2
\end{aligned}
$$

## Many particle, multi level

TBA?
<!-- 
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
## Questions
- [ ] Do we expect the quantum beats if we project the polarization of the light so that we make different polarization of the light indistinguishible?


## Reference
1.  D. A. Steck, Quantum and Atom Optics (2022).
2.  A. Asenjo-Garcia, H. J. Kimble, D. E. Chang, Proc. Natl. Acad. Sci. U.S.A. 116, 25503–25511 (2019).
3.  B. Zhu, J. Cooper, J. Ye, A. M. Rey, Phys. Rev. A. 94, 023612 (2016).
4.  [`CollectiveSpins.jl`](https://qojulia.github.io/CollectiveSpins.jl/dev/descriptions/) theoretical description page.

# Activating this project

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
