# Superfiuidity and phase separation in helium films

Authors: A. N. Berker and David R. Nelson

## Abstract

A vectorial generalization of the Blume-Emery-Griffiths model is proposed to describe superfluidity in films of $^3He-^4He$ mixtures, and is solved by an approximate renormalization scheme due to Migdal. In contrast to bulk mixtures, the line of superfluid transitions is connected to the phase-separation curve by a critical end point. The universal jump of the superfluid density, present in the pure $^4He$ system, is preserved with increasing $^3He$ concentrations $x$ until the critical end point occurs at $x\leq0.12$. At smaller $x$, phase separation causes a kink in the superfluid density versus temperature curve. No tricritical point occurs for any value of the model parameters, although an effectively tricritical phase diagram is obtained in a certain limit. Lines of constant superfluid density bunch up near the effective tricritical point, as predicted by tricritical scaling theory. This treatment also describes superfluidity in pure $^4He$ films in the presence of two-dimensional liquid-gas phase separation. In addition ee calculate the specific heat of the pure $^4He$ system, using the recursion relations of Kosterlitz. This specific heat has a broad maximum above the superfluid transition temperature, corresponding to a gradual dissociation of vortex pairs with increasing temperature.

## 1. INTRODUCTION AND MODEL HAMILTONIAN

In 1971, Blume, Emery, and GriSths introduced a spin-1 Ising model to simulate the thermodynamics of $^3He-^4He$ mixtures along the $\lambda$ line and at the tricritical point. Mean-field treatments of this model display many features of bulk $^3He-^4He$ mixtures, despite the neglect of the continuous rotational symmetry of the degrees of freedom associated with superfluidity. Subsequent work by Riedel and Wegner showed that mean-field theory is in fact correct for tricritical points in all dimensions greater than three. In three dimensions, moreover, the only dependence on the symmetry of the superfluid degrees of freedom is in weak logarithmic corrections to the mean-field predictions. Thus, the success of the Blume-Emery-Grifiths (BEG) model, at least near the tricritical point, is not surprising. Series expansions seem to confirm this general picture.

The BEG model is, however, less appropriate for films of $^3He-^4He$ mixtures. Rigorous mathematics rules out conventional superfluidity, namely, a nonzero order  parameter, in two-dimensional systems. The "superfluid order parameter" of the BEG model is nevertheless expected to be nonzero even in two dimensions, because of the discrete Ising symmetry. Indeed, Monte Carlo and renormalization-group studies of the two-dimensional BEG model yield conventional ordering, although with nonmean-field tricritical exponents.

In this paper, we consider a simple generalization of the BEG model, designed to take into account the continuous rotational symmetry of the superfluid degrees of freedom, and to provide insight into the behavior of films of $^3He-^4He$ mixures. 

Fig.1: tricritical phase diagram of bulk $^3He-^4He$ mixtures from Ref. 32b. The zero-temperature miscibility (coexistence region starts at $x=0.06$ instead of $x=0$) is a quantum-mechanical effect.

BEG used the model Hamiltonian

$$\begin{align}
-\mathcal{H}^{BEG} = J \sum_{\langle i,j \rangle} {s_i s_j} + K \sum_{\langle i,j \rangle} s_i^2 s_j^2 - \Delta \sum_i s_i^2, \tag{1.1}
\end{align}$$

where the spin $s_i$ is located at the site $i$ of some regular lattice, and assumes the discrete values of $0,\pm1$. The first two sums are over nearest-neighbor pairs $\langle i,j \rangle$, $k_B$ and $T$ are the Boltzmann constant and the absolute temperature, and the factor $-1/k_BT$ has been absorbed into the coupling constants $J, K$, and $\Delta$. This can be regarded as a lattice-gas model for mixtures, where the states $s_i=0$ represents the occupation of site $i$ by a $^3He$ atom, and the states $s_i=\pm1$ give the superfluid degrees of freedom when site $i$ is occup,ed by a $^4He$ atom. The generalization of Eq. (1.1) which we have studied in two dimensions in

$$\begin{align}
-\frac{\mathcal{H}}{k_BT} = J \sum_{\langle i,j \rangle} {\vec{s_i} \cdot \vec{s_j}} + K \sum_{\langle i,j \rangle} |\vec{s_i}|^2 |\vec{s_j}|^2 - \Delta \sum_i |\vec{s_i}|^2, \tag{1.2a}
\end{align}$$

where $\vec{s_i}$ is a two-component vector of length unity or zero, i.e.,

$$\begin{align}
\vec{s_i} = (t_i \cos{\theta_i}, t_i \sin{\theta_i}), \tag{1.2b}
\end{align}$$

with $t_i=0,1$ and $0 \leq\ theta_i < 2\pi$. Again $\vec{s_i}=0$ represents occupation by a $^3He$ atom. When site $i$ is occupied by a $^4He$ atom, $\vec{s_i} = (\cos{\theta_i}, \sin{\theta_i})$ reflects the superfluid degrees of freedom, with the appropriate continuous rotational symmetry. Our partition function is

$$\begin{align}
Z = \prod_i{\left( \sum_{t_i=0}^1 \int_0^{2\pi}{\frac{d\theta_i}{2\pi}} \right)e^{-\mathcal{H}/k_BT}}. \tag{1.3}
\end{align}$$

The interpretation of the parameters in (1.2b) is as in the BEG model. The bilinear coupling $k_BTJ$ is a potential promoting superfluid ordering. It is related to bare, areal superfluid density $\rho_0(T)$ by

$$\begin{align}
k_BTJ = \hbar \rho_0(T)/m^2, \tag{1.4}
\end{align}$$

where $m$ is the mass of a $^4He$ atom. The biquadratic coupling $k_bTK$ derives from an isotope effect,

$$\begin{align}
K = K_{33} + K_{44} -2K_{34}, \tag{1.5}
\end{align}$$

when each nearest-neighbor $^\alpha He-^\beta He$ pair is taken to have a classical potential energy $k_BTK_{\alpha\beta}$. The on-site interaction A $k_BT\Delta$ is essentially the difference between the chemical potentials $\mu_3$ and $\mu_4$ of $^3He$ and $^4He$, i.e.,

$$\begin{align}
k_BT\Delta = (\mu_3 - \mu_4) + qk_BT(K_{33} - K_{44}), \tag{1.6}
\end{align}$$

where $q$ is the number of nearest neighbors of a given lattice site. The above is also a model for superfluidity in pure $^4He$ films in the presence of two-dimensional liquid-gas phase separation. Thus, although the remainder of this article employs the terminology $^3He-^4He$ mixtures, the simple reinterpretation as "vacant sites" of "$^3He$ occupied sites" immediately yields a description of superfluidity and condensation in the pure $^4He$ films. The quantity $\langle S_i^2 \rangle$ can then be interpreted as a normalized density of $^4He$ atoms.

In the very negative $\Delta$ limit $(\Delta\rightarrow - \infty)$, the states with any $\vec{s_i}=0$ become negligible in the partition function (1.3), and the planar or XY model is obtained. This can be regarded as a model for superfluidity in pure $^4He$ films, with $\theta_i$ representing the phase of the superfluid order parameter. Fluctuations in the magnitude of this order parameter are not taken into account.

Studies of the two-dimensional XY model have a controversial history. Wegner found from a spin-wave calculation that, at low temperatures, correlations decay algebraically rather than exponentially. High-temperature series expansions by Stanley pointed to the existence of a phase transition at finite temperature, notwithstanding Mermin and Wagner's rigorous proof of the absence of long-range order. A variety of theoretical ideas were subsequently advanced by Berezinskii, Kosterlitz and Thouless, Zittartz, and Luther and Scalapino.

It has become increasingly clear that a simple picture of superfluidity in pure $^4He$ films, due to Kosterlitz and Thouless and to Berezinskii, is substantially correct. One imagines a superfluid state in which long-wavelength phase fluctuations coexist with a dilute gas of bound vortex-antivortex pairs. Although the phase fluctuations prevent true long-range order, correlations decay algebraically, in contrast to the exponential decay in the normal fluid. A transition out of this superfluid state is driven by the dissociation of a finite fraction of the vortex pairs. Kosterlitz, and subsequently Jose et al., have carried out calculations on the two-dimensional planar model which confirm this picture. One striking consequence of these theories is that, at the transition temperature $T_s$ the superfluid density $\rho_s(T)$ in a pure $^4He$ film jumps discontinuously to zero in a universal way. Specifically, it is found that one has

$$\begin{align}
\lim_{T\rightarrow T^-_s}{\frac{\rho_s(T)}{T}} = \frac{2m^2k_B}{\pi\hbar^2} = 3.491x10^{-9} \frac{g}{cm^2.^oK}, \tag{1.7}
\end{align}$$

regardless of the film thickness, substrate, etc. This prediction has been confirmed experimentally in third-sound measurements by Rudnick and collaborators, and in Andronikoshvilli type measurements by Bishop and Reppy.

Here, we study how a systematic dilution with $^3He$ atoms affects the results sketched above. We use an approximate renormalization scheme due to Migdal, which comes remarkably close to the correct result for the pure system. The types of phase diagrams of the dilute system can be anticipated as syntheses of two phenomena: (i) In the pure system, the vortex-unbinding superfluid transition occurs at a temperature $J^{-1} = J^{-1}_s(x=0)$,

$$\begin{align}
J^{-1}_s(x=0) \simeq 1, \tag{1.8}
\end{align}$$

where

$$\begin{align}
x = 1 - \langle s_i^2 \rangle, \tag{1.9}
\end{align}$$

is the $^3He$ concentration. Dilution weakens the effective coupling between the superfluid degrees of freedom, so the transition should occur at lower temperature, i.e., $J_s^{-1}(x)$ decreases as $x$ is increased from 0. (ii) As the $^3He$ concentration $x$ is increased, one expects a phase separation (into two phases, one rich in $^3He$, the other in $^4He$) at temperatures below a criticl temparature 

$$\begin{align}
(J+K)^{-1}_c \sim 1, \tag{1.10}
\end{align}$$

Clearly, the teperature-concentration phase diagrams of our model will strongly depend on $K/J$. 

For $K/J \gg 1$, the phase seperation curve towers over the superfluid transition at x=0. The only readily conceivable phase diagram is shown in Fig. 2. The superfluid transition temperature decreases with increasing x, until it joins the coexistence curve at a critical end point. No multicritical phenomena (such as new exponents) occur at an end point, which is just the point where a higher-order transition gets pre-empted by a first-order line. Such behavior could be expected in, for example, $^3He-Ne$ mixtures, since the "isotope effect" should be rather large in that case. 

For $K/J \ll 1$ on the other hand, another possibility is shown in Fig. 1. The line of superfluid transitions drops down to join the coexistence curve at its tip. This tip is then a tricritical point, which has its own distinctive set of exponents, as evidenced, for example, by the different shapes of the coexistence curves in Figs. 1 and 2. This is the case for bulk $^3He-^4He$ mixtures (to which Fig. 1 applies), whch shows that the isotope effect is weak compared with the superfluid coupling. 

Both types of phase diagrams have been found in the mean-field and renormalization-group investigations of the BEG model. The Migdal renormalization scheme which we use has, in fact, yielded both types of phase diagrams for another two-dimensional system. 

In a given physical situation, as temperature is increased, the Hamiltonian parameters $J$ snd $K$ decrease, because an inverse temperature was absorbed into them [Eq. (1.2a)]. But their ratio $K/J$ is mainly unchanged, and determined by film thickness, susbstrate, etc. This ratio may vary significantly with film thickness. As the thickness is reduced from large values to one or two layers, fluctuations should be more effective in reducing the pure superfluid transition temperature than in diminishing the Ising phase separation $T_c$. Since these characteristic temperatures are determined by Eqs. (1.8) and (1.10) in our model, $K/J$ should increase with decreasing film thickness. It would be rather large for thin films.

Fig. 2: End-point phase diagram of the vectorialized BEG model at $K/J=1$. The dark, full lines indicate higher-order superfluid phase transitions. In (b), constant $\rho_s(T)/k_BT$ curves are drawn with light, full lines: Curve $n$ corresponds to $\rho_s(T)/k_BT = 16 m^2/n\pi\hbar^2$.

The results of our renormalization-group treatment of the vectorialized BEG model (1.2) in two dimensions (on a triangular lattice,which is the close-packing lattice of two dimensions) are summarized as follows. In the space of $J^{-1}, K/J$, and $\Delta$, a first-order surface of concentration discontinuities separates $^3He$-rich and $^4He$-rich phases. This surface terminates in a line of Ising-type critical points.The $^4He$-rich portion of the phase diagram is further divided into superfluid and normal fluid phases by a surface of superfluid phase transitions. The universal jump in the superfluid density, present in the pure system, is preserved at this surface. Our renormalization procedure yields the correlation-length critical exponent $\nu=\infty$ for these superfluid transitions, in agreement with previous theory. The surface of superfluid transitions terminates in a line of critical end points on the first-order surface.

Our temperature-concentration phase diagrams, for large $K/J$, are distinctly of the end-point type, as shown in Fig.2 for $K/J=1$. As $K/J$ is decreased, the end point slides up the coexistence curve, but never actually reaches the phase-separation critical point at the tip. Thus, there is never a tricritical point. The shape of the coexistence curve is always determined by the Ising critical exponent $\beta=1/8$, rather than by some new tricritical phase-separation exponent. As $T$ approaches the phase-separation critical temperature $T$, from below, the coexisting $^3He$ concentrations $x_\pm(T)$ merge as

$$\begin{align}
x_+(T) - x_-(T) \sim (T_c - T)^{1/8}, \tag{1.11}
\end{align}$$

This very flat behavior should be contrasted with that observed in bulk mixtures, wherex $x_{\pm}(T)$  approach each other linearly in $T$. As $K/J$ tends to zero, however, the end point comes very close to the phase-separation critical point. For example, Fig. 3 shows $K/J=0$, the Blume-Capel limit. The end point occurs at $J^{-1}=0.6581$, whereas the phase-separation critical point is at $J^{-1}=0.6615$. The overall phase diagram in Fig. 3 does resemble a tricritical phase diagram. In this case, lines of constant $\rho_s(T)/k_BT$ bunch up as they approach the "effective tricritical point" (formed by the closely situated end point and critical point), as predicted by tricritical scaling theory (Sec. II B).

Fig. 3: Effectively tricritical phased iagram of the vectorialized BEG model at the Blume-Capel plane ($K=0$). The dark, full lines indicate higher-order superfluid phase transitions. In (b), constant $\rho_s(T)/k_BT$ curves are drawn with light, full lines: Curve $n$ corresponds to $\rho_s(T)/k_BT=16m^2/n\pi\hbar^2$.

For any given fixed $K/J$, we find that the line of superfluid phase transitions of the type encountered in pure $^4He$ films, namely, the transition line $T_s(x)$ stretching from $x=0$, is rather short. The first-order transition pre-empts it at 12% of $^3He$ or less, at the critical end point. The $^3He$ concentration at the end point is the maximum amount of $^3He$ which can be included in to a superfluid domain. This limit concentration decreases with increasing $K/J$, being equal to 0.12 for $K/J=0$, and equal to 0.02 for $K/J=1$. These values are quite different from the limit concentration $x=0.67$ of bulk mixtures, which is that of the tricritical point. This difference can be understood qualitatively, by noting that a lesser amount of impurity is needed to drastically weaken the connectivity of interacting superfluid degrees of freedom, when these are arrayed in a lower-dimensional space.

For higher $^3He$ concentrations, the onset of superfluidity in films becomes *nonuniversal*, because of phase separation. As this concentration $x$ is increased from its end-point value, the transition discontinuity in $\rho_s(T)/k_BT$ decreases, but is still finite. This happens in the density interval ($0.12 < x < 0.27$ for $K/J=0$) where phase separation occurs at temperatures higher than the end-point temperature, so that phase separation is initially into a $^3He$-rich phase and a $^4He$-rich normal fluid. The latter undergoes the superfluid transition when the system is further cooled to the end-point temperature. At higher $^3He$ concentrations, the phase-separation temperature is below the end-point temperature, so that the $^4He$-rich phase is always superfluid. In this range, $\rho_s(T)/k_BT$ increases linearly from zero as superfluidity appears at phase separation. Representative superfluid density versus temperature curves are given in Fig. 4. At low $^3He$ concentrations, where the film undergoes the universal superfluid transition, phase separation occurs at a lower temperature. This causes a kink in the superfluid density curve (point $P$ in Fig. 4).

Fig. 4: Superfluid density fraction $\rho_s/\rho_0$ as a function of temperature $(J^{-1})$ at various $^3He$ concentrations $x$, for $K/J=0$. The four types of curve (Sec. III) are represented. Dark dashes indicate that phase separation has occurred. This introduces the kink at point 4P$. Discontinuities at the superfluid transition are shown with light dashes. The locus of minimum superfluidity consists of the universally sloped line segment $U_0 U_E$ the vertical segment at the end-point temperature $J_E^{-1}$, and the $\rho_s/\rho_0=0$ axis between $J_E^{-1}$ and zero temperature.

Fig. 5: Vortex contribution to the specific heat of a pure XY model in two dimensions. The coupling $J$ is the ratio of the nearest-neighbor interaction energy to $—k_BT$. The circle marks the phase transition temperature, where vortex unbinding first occurs. The specific heat has an unobservable essential singularity at this point.

Specific-heat measurements at fixed concentration would reveal any phase separation by a steplike discontinuity. However, such measurements are not a good probe for superfluidity in films. Indeed, one expects only an unobservably weak essential singularity in the specific heat across the superfluid transition. This question is further explored in the Appendix, where the vortex part of the specific heat of the pure XY model is determined. As seen in Fig. 5, this specific heat rises to a maximum at a temperature higher than the superfluid transition temperature $T_s$ with no detectable singularity at $T_s$. The maximum is caused by the gradual dissociation of vortex pairs, beginning at $T_s$. Pairs separated by shorter and shorter distances become unbound with increasing temperatures, until the average separation of bound pairs becomes of the order of magnitude of the vortex core size. The precise shape, height, and position of the maximum is nonuniversal.

Our present model is applicable to films with several atomic layers, provided one deals with quantities averaged over the film thickness. This averaging is especially simple when the correlation length exceeds the film thickness. However, calculations reported here do not take into account any possibly important effect of a concentration gradient (or even phase separation) perpendicular to the film. The Van der Waals attraction to the substrate is the same for $^3He$ and $^4He$ atoms, but the smaller mass and larger zero-point motion of the $^3He$ atoms should cause them to move preferentially to the surface of a thick film. If experiments show that this is of consequence, our model could be further developed by considering 2 two-dimensional lattices stacked on top of each other, with different chemical potentials (i.e., different on-site interactions) at each layer. In thin-film experiments, it is believed that before any superfluidity appears, about one layer of atoms "solidifies" on the substrate. Superfluidity is attributed to subsequent layers. Thus, we can reasonably hope that, in any case, our present calculation is applicable when about one layer is added beyond the solidification layer. It is hard to imagine important vertical differentiation in an effectively single atomic-layer fluid.

The remainder of this paper is organized in the following way: In Sec. II, the renormalization procedure is given. It is argued that $\rho_s(T)/k_BT$ is invariant under this renormalization. Then, starting with the special limits renormalization. Then, starting with the special limits of the pure XY model, the asymptotic first-order region, and spin-(1/2) Ising model, the global Hamiltonian flows are described. In Sec. III, the resulting phase diagrams and, at fixed $^3He$ concentrations, the temperature dependence of the superfluid density are discussed. The details of the specific-heat calculations are relegated to the Appendix. 


## 2. RENORMALIZATION PROCEDURE AND HAMILTONIAN FLOWS

A remarkably simple and powerful method for obtaining recursion relations for comlicated systems in low dimensions has been devised by Migdal. In this section we apply Migdal's method, in a form due to Kadanoff.

### A. Recursion Relations

Consider a slightly more general form of the Hamiltonian (1.2a), namely, 

$$\begin{align}
-\frac{\mathcal{H}}{k_BT} = \sum_{\langle i,j \rangle} {t_i t_j V(\theta_i - \theta_j)} + K \sum_{\langle i,j \rangle} t_i t_j - \Delta \sum_i t_i, \tag{2.1}
\end{align}$$

where $V(\theta)$ is a periodic function with period $2\pi$. Although upon setting

$$\begin{align}
V(\theta) = J \cos{\theta} \tag{2.2}
\end{align}$$

we recover Eq. (1.2a) the specific form (2.2) is not conserved by the renormalization procedure. Therefore, we have to construnct recursion relations for the more general form (2.1). Equation (2.2) is used as initial condition, from which originate the Hamiltonian flows induced by renormalization, in the larger parameter space of Eq. (2.1). (Similarly, if one were solely interested in a system with no biquadratic couplng, $K=0$, the renormalization procedure would generate nonzero $K$, which would then have to be included into the analysis.)

Fig. 6: Migdal transformation for the triangular lattice (Sec. 2 A).

Any renormalization procedure is based on the elimination of a subset of the degresss of freedom inside the partition sum. With our partition sum (1.3), this could be quite a problem. Progress is made, however, by first moving some of the interaction bonds as shown in Fig. 6. One then integrates without difficulty over the isolated one-dimensional degrees of freedom, obtaining effective interactions between the remaining degrees of freedom. The justification of the bond-moving approximation was discussed by Kadanoff. We recall here that this approximation obeys a lower-bound variational pinciple for the free energy. It is known to be fairly accurate at low critcal dimensionalities, which is the case for the $XY$ degrees of freedom in Eq. (2.1). Figure 6 shows the adaptation of bond moving to the triangular lattice, previously used in the study of epitaxial ordering in physisorbed films. 

The bond-moving prescription is completed by considering the treatment of the on-site interactions. One has to decide what fraction of these gets moved with the nearby bonds. For this, a scheme introduced by Emery and Swendsen is employed. The Hamltonian (2.1) is rewritten

$$\begin{align}
-\frac{\mathcal{H}}{k_BT} = \sum_{\langle i,j \rangle} [{t_i t_j V(\theta_i - \theta_j)} - V(0)] * \sum_{\langle i,j \rangle} (t_i - t_j)^2 - [\Delta - \frac{1}{2}q(K+V(0))] \sum_i t_i, \tag{2.3}
\end{align}$$

The first two terms are regarded as bonds, to be moved in their entirety. The last term is regarded as truly on-site, and is not moved. The effect of bond moving should be small at high temperatures, because the moved entities go to zero with $V(\theta)$, $K \sim 1/(k_BT) \rightarrow T$. At low temperatures, the effect of this bond moving should again be small, this time because the configurations with $t_i=t_j$, $\theta_i \simeq \theta_j$ dominate, and the local operators inside the moved entities vanish for these configurations. Finally, note that when we have $V(\theta) = K = 0$, nothing is moved. This means that, along the entire $\Delta$ axis, the partition function is evaluated exactly. This axis acts as an 'anchor' to the approximate evaluations in full $J$, $K$, and $\Delta$ space. These aspects of the Emery-Swendsen scheme are certainly to mitigate the approximation inherent in bond moving. 

The procedure described above has been carried out for a triangular lattice. Starting with a Hamiltonian of the form (2.1), we obtain a new Hamiltonian of the same form, coupling the thinned-out degrees of freedom. In terms of 

$$\begin{align}
u(\theta) \equiv e^{V(\theta)-V(0)}, w \equiv e^{-K-V(0)}, z \equiv e^\Delta \tag{2.4}
\end{align}$$

the parameters of the new Hamiltonian (primed) are given by the simple recursion relations

$$\begin{align}
u'(\theta) = \frac{w^5 z + \int_0^{2\pi}{\frac{d\phi}{2\pi}u^2(\phi)u^2(\theta-\phi)}}{w^5 z + A_4} \tag{2.5a}
\end{align}$$

$$\begin{align}
w' = \frac{(w^3 z + A_2)^2}{(w z + 1)(w^5 z + A_4)} \tag{2.5b}
\end{align}$$

$$\begin{align}
z' = w^9 z \frac{(w z + 1)^6} {(w^3 z + A_2)^6} \tag{2.5c}
\end{align}$$

where

$$\begin{align}
A_n \equiv \int_0^{2\pi} {\frac{d\phi}{2\pi}} u^n(\phi) \tag{2.5d}
\end{align}$$

It is computationally convenient to follow the recursion of the Fourier components $F$ of $u(\theta)$, i.e.,

$$\begin{align}
F\{u(\theta)\} \equiv f(s) = \int_0^{2\pi} \frac{d\theta}{2\pi} e^{i s \theta} u(\theta),
u(\theta) = \sum_{s=-\infty}^{\infty} s^{i s \theta} f(s) .   \tag{2.6}
\end{align}$$

Sums instead of integrals are to be evaluated when 

$$\begin{align}
f'(s \neq 0) = g^2(s) / (w^5 z + A_4), \tag{2.7a}
\end{align}$$

$$\begin{align}
g(s) = F\{u^2 (\theta)\} = \sum_{p=-\infty}^{\infty} f(p) f(s-p), \tag{2.7b}
\end{align}$$

are used instead of Eq.(2.5a). The $s=0$ component can be evaluated from the normalization condition

$$\begin{align}
1 = u'(\theta = 0) = \sum_{s=-\infty}^{\infty} f'(s), \tag{2.8}
\end{align}$$

Also, $f(-s) = f(s)$ is conserved, and Eq. (2.5d) is evaluated as

$$\begin{align}
A_2 = g(s=0), A_4 = \sum_{s=-\infty}^{\infty} g^2(s), \tag{2.9}
\end{align}$$

For the present problem, it was sufticient to keep the $|s| \leq 10$ components, setting to zero the higher ones (which were of negligible magnitude) in order to truncate the infinite sums. However, we did perform spot checks with up to $|s| \leq 40$ nonzero.

### B. Invariance of the superfluid density

In 1996, Josephson argued that the superfluid density in bulk helium near the lambda point would behave as

$$\begin{align}
\rho_s(T) \sim (T_\lambda - T)^\nu, \tag{2.10}
\end{align}$$

where $\nu \simeq 2/3$ is the correlation-length exponent. In $d$ dimensions, this results (which is consequence of the rotational Invariance of the superfluid) regardless

$$\begin{align}
\rho_s(\{K_\alpha\}) = e^{(2-d)l} \rho_s(\{K_\alpha(l)\}), \tag{2.11}
\end{align}$$

where $e^l = b$ is the Kadanoff block size and $\{K_\alpha(l)\}$ are the rescaled parameters. Although the derivation of Eq. (2.12) relies on the existence of long-range order, it does suggest that $\rho_s$ is invariant under a renormalization transformation in precisely $d=2$. 

To explore this question with the vectorialized BEG model, let us introduce a superfluid velocity field 

$$\begin{align}
\vec{v}_{ij} = (\hbar/(m a))(\theta_i - \theta_j) \vec{\delta}_{ij}, \tag{2.13}
\end{align}$$

where $i$ and $j$ are two nearest-neighbor sites, $a$ is the distance between them, and $m$ is the mass of a $4^He$ atom. The velocity (2.13) points along the unit bond-vector $\vec{\delta_{ij}}$. A two-dimensional superfluid density can now be defined in terms of the correlations of the velocity field,

$$\begin{align}
\vec{k_BT/\rho_s} = \frac{2 a^2}{q N} <\sum_{<ij>} \vec{v}_{ij} \cdot \sum_{<kl>} \vec{v}_{kl}>, \tag{2.14}
\end{align}$$

where each sum is over all nearest-neighbor bonds, $(1/2)q N$ is the number of such bonds, and $N$ is the number of sites.

As discussed by Kadanoff, ambiguities arise when Migdal's decimation is used to evaluate crrelations such as Eq. (2.14). A heuristic argument can, however, be constructed for the invariance of $\rho_s/(k_B T)$ under the Migdal procedure i two dimension. Consider a transformation which first 'moves' the bond velocities in Eq. (2.14) together with the bond interactions. For the generalization of the transformation shown in Fig. 6 to block size $b$, one obtains velocities in Eq. (2.14) which are $b$ times the starting Eq. (2.13), i.e.,

$$\begin{align}
\vec{v}_{ij} \rightarrow b \vec{v}_{ij} , \tag{2.15}
\end{align}$$

The terms in each velocity sum of Eq. (2.14) are now simply grouped into subsums $b \vec{v}_{ij} + b \vec{v}_{jk} + b \vec{v}_{pr}$ along lines connecting two sites ($i$ and $r$) which will remain after the renormalization. These linear subsums are immediately performed. For example, for the sites labeled 1, 2, and 3 in Fig. 6, this proceeds as 

$$\begin{align}
b \vec{v}_{1 2} + b \vec{v}_{2 3} = b \frac{\hbar}{m a} (\theta_1 - \theta_3) \vec{\delta}_{1 3} \equiv b^2 \vec{v}'_{1 3}, \tag{2.16}
\end{align}$$

where we recall that the transformed nearest-neighbor distance is $a' = b a$. Dedecorating all the isolated sites such as site 2, Eq. (2.14) becomes

$$\begin{align}
\vec{k_BT/\rho_s} = \frac{2 (b a)^2}{q (b^{-2} N)} <\sum_{<ij>} \vec{v}_{ij} \cdot \sum_{<kl>} \vec{v}_{kl}>', \tag{2.17}
\end{align}$$

where each sum is over all nearest-neighbor bonds of the transformed system, and the thermal average is evaluated using the transformed Hamiltonian. Since we have $a' = b a$ and $N' = b^{-2} N$, the right-hand side of Eq. (2.17) is just $k_B T' / \rho_s'$. Thus, the lattice analogue of Eq. (2.17), for our generalized BEG model in two dimensions, is obtained. Writing $\rho_s / (k_B T)$ as a functional of $V(\theta)$, $K$, and $\Delta$, we have

$$\begin{align}
\frac{\rho_s}{k_BT}\{V(\theta),K,\Delta\} = \frac{\rho_s}{k_BT} \{V'(\theta),K',\Delta'\}, \tag{2.18}
\end{align}$$

The superfluid density, normalized by absolute temperature, is indeed left  unchanged by the Migdal transformation in two dimensions.

According to the above argument, renormalization-group trajectories are also lines of constant $\rho_s/k_BT$. This invariance property can be used to obtain superfluid densities, provided that Hamiltonians eventually interate into regions where Eq. (2.14) is easily evaluated. As discussed in Sec. 2 E, in the present calculation, Hamiltonians of the superfluid phase iterate onto the Villain model line of fixed points. Equation (2.14) is then readily evaluated by extending the angle integrations to $\pm \infty$, yielding

$$\begin{align}
(\rho_s/k_BT) \simeq (m^2/\hbar^2)J_\nu, \tag{2.19}
\end{align}$$

where $J_\nu$ parametrizes the Villain model as in Eq. (2.25). This result, used in conjunction with Eq. (2.18), produced the lines of constant $\rho_s/k_BT$ in the phase diagrams of Figs. 2,3,7, and 8, and the superfluid density curves of Fig. 4. 

The invariance of $\rho_s/k_BT$ would have interesting consequences near a tricritical point, should one actually exist. In that case, $\rho_s/k_BT$ would become a function of just two relevant thermodynamic fields, say $K_1$ and $K_2$. The phenomenological and renormalization-group formulations of tricritical scaling would then assert thet

$$\begin{align}
\frac{\rho_s}{k_BT}(K_1,K_2) = \frac{\rho_s}{k_BT}(b^{\lambda_1}K_1, b^{\lambda_2}K_2), \tag{2.20}
\end{align}$$

where $K_1=K_2=0$ corresponds to the critical point. 

Fig. 7: End-point phase diagram of the vectorialized BEG model at $K/J=1$. The dark, full line and the line of open circles (ooo) respectively indicate higher-order and first-order transitions. The dark circle (o) marks the isolated critical point, and $E$ labels the critical end point. Constant $\rho_s(T)/k_BT$ curves are drawn with light, full lines: Curve $n$ corresponds to $\rho_s(T)/k_BT = 16m^2 / n\pi\hbar^2$.

This implies that $\rho_s(T)/k_BT$ is a function of only $K_1/K_2^\phi$, with $\phi=\lambda_2/\lambda_1$, i.e.,

$$\begin{align}
\frac{\rho_s}{k_BT}(K_1,K_2) = R \frac{K_1}{K_2^\phi}, \tag{2.21}
\end{align}$$

We are led to the conclusion that the finite limiting value of $\rho_s/k_BT$ at the tricritical point is indeterminate, depending instead on the path of approach! Lines of constant $\rho_s/k_BT$ converge onto the tricritical point. Although we found no tricritical fixed point in our calculation, this convergence is exhibite at the scale of Fig.8(a). The constant $\rho_s/k_BT$ curves labeled 8 through 5 apparently meet at point $T_e$ and have $\rho_s/k_BT$ values of $2m^2/\pi\hbar^2$ through (3.2) $m^2/\pi\hbar^2$.

Fig. 8: (a) Effectively tricritical phase diagram of the vectorialized BEG model at the Blume-Capel plane $(K=0)$. The dark, full line and the line of open circles (ooo) respectively indicate higher-order and first-order transitions. At the scale of this figure, these lines meet at the effective tricritical point $T_e$. Constant $\rho_s(T)/k_BT$ curves are drawn with light, full lines: Curve $n$ corresponds to $\rho_s(T)/k_BT = 16m^2/n\pi\hbar^2$. As predicted by tricritical scaling theory (Sec. 2 B), curves 5-8 apparently converge onto the effective tricritical point. (b) Detailed view of the region $T_e$. [In going from (a) to (b), the horizontal and vertical scales were respectively blown up by factors of 333 and 125, so that (b) actually covers less area than the black point at $T_e$ in (a).] The end-point structure is seen here: (o) nd (*) respectively mark the critical point terminating the first-order line and the end point terminating the higher-order line.

### C. Pure XY limit

For $\Delta \rightarrow -\infty$, Eq. (2.5a) we recover Migdal's recursion relation for the generalized XY Hamiltonian

$$\begin{align}
\frac{\mathcal{H}^{XY}} {k_BT} = \sum_0^{2\pi} {[V(\theta_i - \theta_j) - V(0)]} , \tag{2.22}
\end{align}$$

namely,


$$\begin{align}
u'(\theta) = A_4^{-1} \int_0^{2\pi} \frac{d\phi}{2\pi} {u^2(\phi) u^2(\theta - \phi)} , \tag{2.23}
\end{align}$$

or in terns of Fourier components


$$\begin{align}
f'(s) = A_4^{-1} \left[ \sum_{p=-\infty}^{\infty} {f(p) f(s-p)} \right]^2 , \tag{2.24}
\end{align}$$

This recursion relation was studied in some detail by Migdal and by Jose et al. As discussed in detail by the latter authors, a low temperature, but otherwise arbitrary $V(\theta)$, relaxes after a few iterations to a Villain potential with Fourier components

$$\begin{align}
f_{\nu}^*(s) = e^{-s^2/2J_{\nu}} \left[ \sum_{p=-\infty}^{\infty} {e^{-p^2/2J_{\nu}}} \right]^{-1} , \tag{2.25}
\end{align}$$

Herer again, the coupling constant $J_{\nu}$ should be thought of as an inverse tempreature. 

At low temperatures, the Villain potential (2.25) is very nearly invariant under the renormalization procedure, as can be checked by substituting Eq. (2.25) into Eq. (2.24), and approximating the sums by integrals. In fact, a "quasifixed line" of functions of the form (2.25) is found numerically for $J_{\nu}^{-1} \leq \frac{1}{2}\pi$. Unfortunately, a small drift toward high temperature is always present, so that all finite-temperature initial conditions eventually reach an infinite-temperature fixed point at $V(\theta) = 0$. However, for example, an initial potential (2.25) with $J_{\nu}^{-1} = \frac{3}{2}(\frac{1}{2}\pi)$ takes 11 iterations to reach $V(\theta) \sim 10^{-38}$, whereas over 12000 iterations are necessary for an initial $J_{\nu}^{-1} = \frac{1}{2}(\frac{1}{2}\pi)$. For still lower temperatures, Eq. (2.25) is numerically indistinguishable from a fixed point. Since more rigorous treatments of the XY model to produce genuine fixed-line behaviour terminating at $J_{\nu}^{-1} \simeq \frac{1}{2}\pi$; the drift here must be regarded as a defect of Migdal's renormalization procedure. In this work, we effectively ignore this small drift toward high temperature, and consider

$$\begin{align}
f(s) = f_{\nu}^*(s) , J_{\nu}^{-1} \leq \frac{1}{2}\pi , \Delta \rightarrow -\infty , \tag{2.26}
\end{align}$$

to be a line segment of fixed points ($S_0 S_1$ in Fig. 9). The initial coupling given in Eq. (2.2), namely, $V(\theta) = J \cos{\theta}$, maps onto this fixed line for $J^{-1} \leq 1.025$.

Finally, in this $\Delta \rightarrow -\infty$, states with any $t_i = 0$ are negligible şn the partition sum (1.3), and the biquadratic coupling K becomes an additive constant in the Hamiltonian. Concurrently, its recursion relation (2.5b) is slaved to $V(\theta)$, i.e.,

$$\begin{align}
w' = A_2^2 / A_4 , \tag{2.27}
\end{align}$$

On the other hand, the coupling $\Delta$ has a crucial role: its renormalization determines the stability of the $\Delta \rightarrow -\infty$ limit. The recursion relation (2.5c) reduces to

$$\begin{align}
e^{\Delta'} = z' = w^9 z A_2^{-6} . \tag{2.28}
\end{align}$$

Fig. 9: Global Hamiltonian flows, schematically shown (Sec.II E) Flows through higher-order superfluid transitions, first-order transitions, and Ising-type critical points are respectively indicated by dark full lines, lines of open circles (ooo), and a line of dark circles (***). The flow drawn with alternating dashes and open circles (o—o) is through critical end points. The flow drawn with alternating dashes, open and dark circles (o—*—o) is where the Ising-type critical points and the end points come very close together to form effective tricritical points. Fixed points and fixed lines are indicated by (^). The paramagnetic fixed line is not shown.

At the fixed line (2.26), we can substitute Eq. (2.27) into Eq. (2.28) to obtain

$$\begin{align}
z'/z = A_2^{12} A_4^{-9} . \tag{2.29}
\end{align}$$

Evaluating $A_n$ for the Villain potential (2.25), approximating sums by integrals we have

$$\begin{align}
 \frac{z'}{z} = \left[ \frac{2}{\pi} \frac{1}{J_{\nu}} \right]^{3/2} . \tag{2.30}
\end{align}$$

Thus, the $\Delta \rightarrow -\infty$ limit is stable at the fixed-line segment (2.26), which is therefore expected to be important in the global phase diagram of the dilute system. From this calculation, we have $\Delta'=\Delta$ at exactly $J_{\nu}^{-1} = \frac{1}{2}\pi$, so the checmical potential becomes a marginal variable, and the pure XY limit reverses stability at exactly the edge of the fixed-line segment (2.26). At $J_{\nu}^{-1} > \frac{1}{2}\pi$, $\Delta$ is renormalized to kess negative values. The above is also precisely the stability behavior of the vortex-core energy $\ln{y}$ which arises in analyti treatments of the pure XY model. Indeen, one would expect $^3He$ atoms to cluster around vortex cores, thus contributing to the core energy. 

### D. First-order fixed line and Ising behavior

In the limit $\frac{1}{2}qK \sim \Delta \rightarrow \infty$, the pure XY recursion relation (2.23) is again recovered, and the preceding discussion applies. A new fixed-line segment ($F_0 F_1$ in Fig. 9) is located at

$$\begin{align}
f(s) = f_{\nu}^*(s), J_{\nu}{-1} \leq \frac{1}{2}\pi , \tag{2.31a}
\end{align}$$

$$\begin{align}
\Delta \rightarrow \infty , \tag{2.31b}
\end{align}$$

$$\begin{align}
\frac{1}{2}qK - \Delta - \ln{A_4} = 0 , \tag{2.31c}
\end{align}$$

where $A_4$ is evaluated by substituting into Eq. (2.9) the Villain potential (2.25), and is approximately $(8\pi J_{\nu})^{1/2}$. This fixed-line segment has one stable direction, $\Delta' \simeq 2\Delta \rightarrow \infty$, and one unstable direction involving Eq. (2.31c). For small deviations from zero we have

$$\begin{align}
\frac{1}{2}qK' - \Delta' - \ln{A_4} = b^d (\frac{1}{2}qK - \Delta - \ln{A_4}) , \tag{2.32}
\end{align}$$

where $b=2$ is the length rescaling factor and $d=2$ is the lattice dimensionality. Equation (2.32) is a necessary condition for a fixed point whose domain of attraction is a locus of first-order phase transitions, as noted by Nauenberg and Nienhuis. Indeed, our density calculations reveal in the present case a locus of discontinuities in concentration.

These first-order phase transitions can be located also by noting, in the $(\frac{1}{2}q)K \sim \Delta \rightarrow \infty$ limit, the crossing of the energies of the two dominant configurations, $\{ \vec{s}_i = \vec{s}_j = ... = \vec{s}_N \neq 0 \}$ and $\{ \vec{s}_i = \vec{s}_j = ... = \vec{s}_N = 0 \}.$. This exercise would immediately yield the coeificient $\frac{1}{2}q$ in Eq. (2.31c). The logarithmic correction $\ln{A_4}$ is a bit more subtle. It derives from the fact that the former configuration (aligned nonzero spins) actually has zero weight in the partition sum, but occurs together with a continuum of spin-wave configurations with infinitesimal energy variations. This logarithmic term does not appear in the original BEG model, which has only discrete excitations. For the present vectorialized model, are normalization procedure which equipartitions the on-site interactions onto the bonds, in contrast to the Emery-Swendsen scheme, gives the same logarithmic term.

As the coupling between the XY degrees of freedom, e.g., $J_{\nu}$ in Eq. (2.25) or $J$ in Eq. (2.2), is weakened, the fixed-line segment (2.31) becomes unstable to the isolated first-order fixed point ($F$ in Fig. 9)

$$\begin{align}
V(\theta) = 0, i.e., f(s) = \delta_{s0} , \tag{2.33a}
\end{align}$$

$$\begin{align}
\Delta \rightarrow \infty , \tag{2.33b}
\end{align}$$

$$\begin{align}
\frac{1}{2}qK - \Delta = 0 . \tag{2.33c}
\end{align}$$

Equation (2.33b) is again stable, $\Delta' \simeq 2\Delta$, and deviations from Eq. (2.33c) again satisfy the Nauenberg-Nienhuis condition. $V(\theta)=0$ is completely stable. This fixed point occurs in the region $V(\theta) where the BEG model reduces to a spin-(1/2) Ising model. Defining a new variable 

$$\begin{align}
\sigma_i \equiv 2t_i -1 , \tag{2.34}
\end{align}$$

the $V(\theta) = 0$ Hamiltonian (2.1) becomes

$$\begin{align}
-\frac{\mathcal{H}'} {k_BT} = \frac{K}{4} \sum_{\langle ij \rangle } {\sigma_i \sigma_j} + \frac{1}{2}(\frac{1}{2}qK-\Delta) \sum_i{\sigma_i} ; \sigma_i = \pm1 . \tag{2.35}
\end{align}$$

The line $(\frac{1}{2}q)K=\Delta$ corresponds to zero magnetic field in this Ising model, and is closed under our renormalization procedure. The large $K$ (i.e., low temperature) segment of this line is a first-order boundary between $\langle t_i \rangle \gtrless \frac{1}{2}$ phases. Initial conditions on this segment flow, under successive renormalizations, to the fixed point (2.33). The small $K$ (high temperature) segment forms part of a smooth trend (no phase transition) across $\langle t_i \rangle = \frac{1}{2}$, and flows to a high-temperature fixed point at $V(\theta) = K = \Delta = 0$. The low- and high-temperature segments are separated by the critical fixed point ($C$ in Fig. 9) at 

$$\begin{align}
V(\theta)=0, K_c^* = 1.2188, \Delta_c^* = (\frac{1}{2}q)K_c^* , \tag{2.36}
\end{align}$$

Within our approximate treatment, unstable (relevant) deviations from this fixed point have the eigenvalue exponents $\lambda_H=1.797$ for the magnetic field direction and $\lambda_T=0.747$ for the thermal direction. These numbers are to be compared with the exact results

$$\begin{align}
K_c=\ln{3}=1.0986, \lambda_H=\frac{15}{8}, \lambda_T=1 .
\end{align}$$

### E. Global Hamiltonian flows

The global Hamiltonian flows, which connect the special limits of Sec. II C and II D, are given schematically in Fig. 9. In this figure, the axis parameter

$$\begin{align}
J_{eff} \equiv -\frac{d^2 V(\theta)}{d \theta^2} |_{\theta=0}, \tag{2.37}
\end{align}$$

introduced in Ref. 22, measures the coupling between the XY degrees of freedom. It equals $J$ for the initial cosine potential (2.2), and it is approximately (sums replaced by integrals) $J_{\nu}$ for the Villain potential (2.25). The choice of the other axis parameter, 

$$\begin{align}
D \equiv (\Delta + \ln{A_4})/K, \tag{2.38}
\end{align}$$

is motivated by the asymptotic first-order conditions (2.31c) and (2.33c).

In Fig. 9, the domain of attraction $C_0 C$ of the Ising-critical fixed-point $C$ [Eq. (2.36)] terminates the first-order domains of the fixed-line segment $F_0 F_1$ [Eq. (2.31)] and of the fixed-point $F$ [Eq. (2.33)]. Within these combined first-order domains, $F_0 F_1 \ (F)$ is the terminus of flows from regions of strong (weak) XY coupling $J_{eff}$.

The fixed-line segment $S_0 S_1$ [Eq. (2.26)] of the pure XY model is the sink for the superfluid region of the dilute system. This superfluid region is at small $J_{eff}^{-1}$, and to the left of the figure. It is bounded by two surfaces: (i) Flows terminating along the whole length of $F_0 F_1$. These flows constitute a boundary of first-order phase transitions. (ii) Flows terminating at $S_1$, the edge of the fixed-line segment. Separating $\Delta \rightarrow -\infty$  limits of opposite stability (Sec. II C), these flows are a renormalization-group separatrix. They constitute a boundary of higher-order phase transitions, of the same type as in the pure XY model. Since $\Delta$ is marginal [Eq. (2.30)] at $S_1$, the correlation-length critical exponent is

$$\begin{align}
\nu = \infty . \tag{2.39}
\end{align}$$

The high-order boundary mentioned above terminates on the first-order domains in a line $E F_1$ of critical end points, between the domains of $F_0 F_1$ and $F$. Although for small $K/J_{eff}$ the line $C_0 C$, there is always a strip in between them which flows to the Ising first-order fixed-point $F$. Unlike the renormalization-group treatments of the discrete BEG model, no tricritical fixed point occurs. 

Finally, not shown in Fig. 9, a "paramagnetic" fixed line at $V(\theta) = K = 0$, $\Delta$ arbitrary, is the sink for the normal-fluid regions.

