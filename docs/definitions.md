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

## Green's function

The imaginary part of the Green's function defined in `src/base.jl` approaces to 2/3 at $\mathbf{r}_{ij} = (0, 0, 0)$. See `test/GreenTensor_test.jl`. To account this normalization, we put 3/2 and 3/4 for calculation J and Γ. See for example, `src/many_atoms_four_level.jl`.

```julia
    # J and Γ matrix
J(i, j, p, q, kl::CartesianIndex) = - 3/4*Γ_kl[kl]*ϵ_q(p)*real.(GreenTensor(positions[i] - positions[j], knorm[kl])) * adjoint(ϵ_q(q))
Γ(i, j, p, q, kl::CartesianIndex) =   3/2*Γ_kl[kl]*ϵ_q(p)*imag.(GreenTensor(positions[i] - positions[j], knorm[kl])) * adjoint(ϵ_q(q))
```



## Reference

1. A. Asenjo-Garcia, H. J. Kimble, D. E. Chang, Proc. Natl. Acad. Sci. U.S.A. 116, 25503–25511 (2019).
