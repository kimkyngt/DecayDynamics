# Definitions

## Atomic rasing operator

We consider magnetic sub-level contributions similar to equation (2) in the reference [1]. However, this expression seems to make total decay rate to be dependant on the angular momentum quantum number of the excited state. See `src/Gamma_0_test.jl` to see its dependence. To handle the total decay rate easily, we add a factor $\sqrt{2F_e + 1}$ to the definition of the atomic rasing operator(equation (2) in [1]).

$$\begin{gathered}
\hat{\Sigma}_q = \sum_{m_g = -F_g}^{F_g}  \sqrt{2F_e + 1} (-1)^{F_g-m_g}\left(\begin{array}{llc}
F_g & 1 & F_e \\
-m_g & q & m_g-q
\end{array}\right) |F_e m_g-q\rangle \langle F_g m_g |  \\
= \sum_{m_g = -F_g}^{F_g} (-1)^{2m_g - q} |F_g m_g \rangle \langle F_g, m_g; 1, -q| F_e, m_g - q \rangle \langle F_e m_g - q |
\end{gathered}$$

The other parts are essentially the same as [1].

$$\begin{gathered}
\mathscr{H}=\hbar \sum_{i, j=1}^N \sum_{q, q^{\prime}=-1}^1 J_{i j q q^{\prime}} \hat{\Sigma}_{i q}^{\dagger} \hat{\Sigma}_{j q^{\prime}}, \\ \mathscr{L}[\rho]=\sum_{i, j=1}^N \sum_{q, q^{\prime}=-1}^1 \frac{\Gamma_{i j q q^{\prime}}}{2}\left(2 \hat{\Sigma}_{j q^{\prime}} \rho \hat{\Sigma}_{i q}^{\dagger}-\hat{\Sigma}_{i q}^{\dagger} \hat{\Sigma}_{j q^{\prime}} \rho -\rho \hat{\Sigma}_{i q}^{\dagger} \hat{\Sigma}_{j q^{\prime}}\right)
\end{gathered}$$

$$
\begin{aligned}
& J_{i j q q^{\prime}}=-\frac{\mu_0 \omega_0^2}{\hbar}|\wp|^2 \hat{\mathbf{e}}_q \cdot \operatorname{Re} \mathbf{G}\left(\mathbf{r}_i, \mathbf{r}_j, \omega_0\right) \cdot \hat{\mathbf{e}}_{q^{\prime}}^* \\
& \Gamma_{i j q q^{\prime}}=\frac{2 \mu_0 \omega_0^2}{\hbar}|\wp|^2 \hat{\mathbf{e}}_q \cdot \operatorname{Im} \mathbf{G}\left(\mathbf{r}_i, \mathbf{r}_j, \omega_0\right) \cdot \hat{\mathbf{e}}_{q^{\prime}}^*
\end{aligned}
$$

## Green's function, $J$ and $\Gamma$

The imaginary part of the Green's function defined in `src/base.jl` approaces to 2/3 at $\mathbf{r}_{ij} = (0, 0, 0)$. See `test/GreenTensor_test.jl`. To account this normalization, we put 3/2 and 3/4 for calculation J and Γ. See for example, `src/many_atoms_four_level.jl`. Here, indices `kl` is for the different transitions.

```julia
    # J and Γ matrix
J(i, j, p, q, kl::CartesianIndex) = - 3/4*Γ_kl[kl]*ϵ_q(p)*real.(GreenTensor(positions[i] - positions[j], knorm[kl])) * adjoint(ϵ_q(q))
Γ(i, j, p, q, kl::CartesianIndex) =   3/2*Γ_kl[kl]*ϵ_q(p)*imag.(GreenTensor(positions[i] - positions[j], knorm[kl])) * adjoint(ϵ_q(q))
```


## Reference

1. A. Asenjo-Garcia, H. J. Kimble, D. E. Chang, Proc. Natl. Acad. Sci. U.S.A. 116, 25503–25511 (2019).
