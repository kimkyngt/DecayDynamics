# Dipole operator with magnetic sublevels

Following [1] equation (7.396), atomic lowering operator for $F_g \rightarrow F_e$ is

$$\begin{gathered}
\hat{\Sigma}_q = \sum^{F_e}_{m_e=-F_e} |F_g m_e + q \rangle \langle F_g m_e + q | F_e m_e ; 1 q \rangle \langle F_e m_e | \\
 = \sum^{F_g}_{m_g=-F_g} |F_g m_g \rangle \langle F_g m_g | F_e m_g - q ; 1 q \rangle \langle F_e m_g - q |.
\end{gathered}$$

The CGC can be expressed as

$$\begin{gathered}
\langle F_e m_g - q ; 1 q |  F_g m_g\rangle = (-1)^{1+q} \sqrt{\frac{2F_g + 1}{2F_e + 1}} \langle 1 - q; F_g m_g | F_e m_g -q\rangle \\
= (-1)^{1+q + 1 +F_g - F_e} \sqrt{\frac{2F_g + 1}{2F_e + 1}} \langle  F_g m_g ;1 - q | F_e m_g -q\rangle \\
= (-1)^{1+q + 2(1 +F_g - F_e)} \sqrt{\frac{2F_g + 1}{2F_e + 1}} \langle  F_g - m_g ;1 + q | F_e -(m_g -q)\rangle
\end{gathered}$$

The CGC can be expreesed in terms of [3-j symbols](https://en.wikipedia.org/wiki/Clebsch–Gordan_coefficients),

$$
\begin{aligned}
\left\langle j_1 m_1 j_2 m_2 \mid J M\right\rangle & =(-1)^{-j_1+j_2-M} \sqrt{2 J+1}\left(\begin{array}{ccc}
j_1 & j_2 & J \\
m_1 & m_2 & -M
\end{array}\right)
\end{aligned}
$$

$$\begin{gathered}
\langle F_e m_g - q ; 1 q |  F_g m_g\rangle = (-1)^{1+q + 2(1 +F_g - F_e)} \sqrt{\frac{2F_g + 1}{2F_e + 1}}  (-1)^{-F_g + 1 - m_g + q } \sqrt{2F_e+1}\left(\begin{array}{ccc}
F_g & 1 & F_e \\
-m_g & q & m_g-q
\end{array}\right)  \\
= (-1)^{-F_g-m_g}\sqrt{2F_g + 1}\left(\begin{array}{ccc}
F_g & 1 & F_e \\
-m_g & q & m_g-q
\end{array}\right)
\end{gathered}$$

Factor $\sqrt{2F_g + 1}$ difference from [2] is likely to be came from the different definition of the reduced dipole matrix element. In [1], the following definition has been used

$$
\langle j m | T^{(k)}_q | j' m'\rangle = \langle j' m' k q | j m \rangle \langle j \| T^{(k)} \| j'\rangle.
$$

Whereas, in [2] an alternative version seems to be used

$$
 \langle j m | T^{(k)}_q | j' m'\rangle
  = \frac{(-1)^{2 k} \langle j' m' k q | j m \rangle \langle j \| T^{(k)} \| j'\rangle_{\mathrm{R}}}{\sqrt{2 j + 1}}
  = (-1)^{j - m}
    \begin{pmatrix}
      j & k &  j' \\
      -m & q & m'
    \end{pmatrix} \langle j \| T^{(k)} \| j'\rangle_{\mathrm{R}}.
$$

We follow the convention from [1].

With $\Sigma_q$, we write atom-field interaction Hamiltonian as

$$
H_{\mathrm{AF}}=\frac{\hbar}{2} \sum_q\left[\Omega_q^* \Sigma_q e^{i \omega t}+\Omega_q \Sigma_q^{\dagger} e^{-i \omega t}\right]
$$

With rabi frequency defined as follows.
$$
\Omega_q= - \frac{2\left\langle J_{\mathrm{g}}\|\mathbf{d}\| J_{\mathrm{e}}\right\rangle E_{q}^{(+)}(0)}{\hbar}
$$

## Reference

1. D. A. Steck, Quantum and Atom Optics (2022).

2. A. Asenjo-Garcia, H. J. Kimble, D. E. Chang, Proc. Natl. Acad. Sci. U.S.A. 116, 25503–25511 (2019).
