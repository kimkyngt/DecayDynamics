# Atomic rasing operator

We consider magnetic sub-level contributions similar to equation (2) in the reference [1]. However, this expression seems to make total decay rate to be dependant on the angular momentum quantum number of the excited state. See `src/Gamma_0_test.jl` to see its dependence. To handle the total decay rate easily, we add a factor $\sqrt{2F_e + 1}$ to the definition of the atomic rasing operator(equation (2) in [1]).

$$\begin{gathered}
\hat{\Sigma}_q = \sum_{m_g = -F_g}^{F_g}  \sqrt{2F_e + 1} (-1)^{F_g-m_g}\left(\begin{array}{llc}
F_g & 1 & F_e \\
-m_g & q & m_g-q
\end{array}\right) |F_e m_g-q\rangle \langle F_g m_g |  \\
= \sum_{m_g = -F_g}^{F_g} (-1)^{2m_g - q} |F_g m_g \rangle \langle F_g, m_g; 1, -q| F_e, m_g - q \rangle \langle F_e m_g - q |
\end{gathered}$$

The other parts are essentially the same as [1].

## Reference

1. A. Asenjo-Garcia, H. J. Kimble, D. E. Chang, Proc. Natl. Acad. Sci. U.S.A. 116, 25503–25511 (2019).
