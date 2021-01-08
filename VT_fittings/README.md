# Fitting of volume _versus_ Temperature and estimation of $\alpha(T)$

A equation was fitted to the volume versus temperature data contained in the .dat files, which, in turn, was obtained from Molecular Dynamics simulations.

The equation used to describe the volume as a function of temperature is a scaled summation of Einstein-like terms,

$$
V(\tau) = a_0 \sum_{i=1}^N \frac{a_i b_i}{e^{b_i/\tau}-1},
$$

where $\tau=\frac{T}{1000}$, $T$ is the temperature and $a_i$, $b_i$ are the parameters of the equation. In all cases, we have truncated the equation to N=3. The parameters are determined using subsequent fitting routines from the [LMFit](https://lmfit.github.io/lmfit-py/index.html) package, which includes the [AMPGO solver](http://infinity77.net/global_optimization/ampgo.html) and Nelder-Mead. The fitted parameters are then used to calcualte the volumetric coefficient of thermal expansion $alpha$ from its thermodynamical definition,

$$
\alpha (T) = \frac{1}{V(T)}\frac{\partial V}{\partial T},
$$

which yealds

$$
\alpha(\tau) = \sum_{i=1}^{N}\frac{a_i b_i e^{b_i/\tau}}{\tau^2\left(e^{b_i/\tau}\right)^2}
$$

All the fittings and its plots were made in Jupyter notebooks contained in the relevant folders. For reference purposes, the plots of $\alpha$ contains data points calculated from finite differences of the volume values.
After each fit a *.fit file is created with columns containing the data and fitted values.

A summary file containing two columns of alpha at 100 K and 300K for each CRN is also generated. This file is used for the generation of the box plots presented in the main paper.

Note: The AMPGO solver is a stochastic algorithm, so one might get slightly different results for an individual fitting. This does not changes the collective behaviour.