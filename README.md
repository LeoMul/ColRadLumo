# ColRadLumo

This code internally uses [ColRadPy](https://github.com/johnson-c/ColRadPy) ([paper citation, Curt Johnson (2019)](https://doi.org/10.1016/j.nme.2019.01.013)) as a library to solve the statistical rate equations. The resulting data is then post-processed into astrophysical absolute emission lines. The code is designed to be ran interactively in a JuPyter notebook. The input data is specified by user in the adf04 format.

**If you use this code**
Please cite Curt's original [code](https://github.com/johnson-c/ColRadPy) and [paper](https://doi.org/10.1016/j.nme.2019.01.013), and this repository if possible.


Credit to ColRadPy: https://github.com/johnson-c/ColRadPy
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 

*Features*
1. Calculation of astrophysical luminosities,
$$
    L_{j \to i} = \frac{hc}{\lambda_{j \to i}}   \frac{n_e\text{PEC}_{j\to i } }{\sum_i N_i} \frac{M_{\text{ion}}}{m_{\text{ion}}} ~~\text{[erg s$^{-1}$]}
$$ where $m_{\text{ion}}$ is the nuclear mass of your atomic species - which is taken from a lookup table in `atomic_masses.py`. This feature requires you to have set the atomic symbol correctly in the adf04 file.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 
2. Spectral synthesis, assuming some FWHM velocity $\beta$ and central wavelength 
$\lambda_0$ (nm) - 
$$
\sigma_\lambda = \lambda_0 \beta /2.355 ~~\text{[nm]} \\
L^{\lambda}_{j \to i} = \frac{1}{10\sigma_\lambda \sqrt{2\pi}} \exp\left( -\frac{1}{2} \left(\frac{\lambda_0 - \lambda}{\sigma_\lambda}\right)^2\right) ~~\text{[erg s$^{-1}$ Å$^{-1}$]}
$$
where the additional factor of 10 changes the units to Angstrom.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 
3. Production of mass-density-temperature contour plots for diagnosis of spectral lines. (public code coming)
4. Isolate strongest emission lines 
5. Basic self-consistent inclusion of Sobolov depth,
$$
\begin{align*}
    \beta_{i\to j} &= \frac{1}{\tau_{i\to j}} \left(1-e^{-\tau_{i\to j}} \right), \text{with}\\
    \tau_{i\to j} &= \frac{A_{i\to j} \lambda^3 n_j t}{8 \pi}\left( \frac{g_i}{g_j} -  \frac{n_i}{n_j}\right)
\end{align*}
$$
where the level densities $n_j \propto N_j$ are calculated either at some velocity and explosion time, or as a fraction of the electron density. At the first iteration, we take $\beta_{i\to j}^0 = 1$ and calculate the populatons $N_j$. This is used to calculate a new estimate to $\beta_{i\to j}^1$ We then gradually add on the opacity at iteration $k$ ala 
$$
\begin{equation*}
    \beta_{i\to j}^k \to \frac{1}{2} \beta_{i\to j}^k + \frac{1}{2} \beta_{i\to j}^{k-1}
\end{equation*}
$$
as is commonly done in self-consistent iterations. Internally, the ColRadPy class has its A-value arrays updated according to $$A_{i \to j} \to A_{i \to j}\beta^k_{i \to j,}$$ until convergence. This presently is only implemeneted for a single density and temperature point - in general the Sobolov escape factors are themselves dependent on temperature and density. 

*Works using this code*:
1. [Mulholland, McNeill et al (2024)](https://doi.org/10.1093/mnras/stae2331)
2. [McCann, Mulholland, Xiong et al (2025)](https://doi.org/10.1093/mnras/staf283).
3. [Mulholland et al (2026)](https://doi.org/10.1093/mnras/stag237)


*Current todo's*
1. Presently - might make module a pre requisite -  change path to your installation of ColRadpy in colradlumo_library.py
2. Test and upload mass-estimate contours
3. Figure out how to self-consistently iterate Sobolov for *all* density and temperature points. 

