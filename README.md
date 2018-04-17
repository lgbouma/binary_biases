The code in this repository supports the results presented by Bouma, Masuda,
and Winn in a paper entitled "Biases in Planet Occurrence Caused by Unresolved
Binaries in Transit Surveys".

## Dependencies

For calculations:
  * `python3`
  * `numpy`
  * `scipy`
  * `pandas`

For plots:
  * `matplotlib`
  * `src/brokenaxes.py` [forked from here](https://github.com/bendichter/brokenaxes)

One easy way to install all of the above is with the [Anaconda python distribution](https://conda.io/docs/user-guide/install/index.html).

## To reproduce the main plots

`src/integrate_for_apparent_rate_density.py`: integrate the general equation
for apparent rate density. The output from this code is written to
`/data/numerics`.

`src/run_integrator.sh`: shell script that runs all the models used in the
paper.

`src/plot_integrated_rate_densities_vs_radius.py` make plots of apparent radius
density as a function of planet radius.

## Other goodies

`src/numerical_models.py`: run Monte Carlo surveys for (A) twin binary,
varying planet populations, and for (B) fixed planet radius, varying binary
mass ratio populations. The routine behind this is explained in more depth in
`monte_carlo_readme.md`, and has only been tested for these two cases.
