This code supports the results presented by Bouma, Masuda, and Winn in a paper
called

BINARITY’S BIASES ON TRANSIT SURVEY OCCURRENCE RATES

## To reproduce the plots:

`src/integrate_for_apparent_rate_density.py` integrates the general equation
for apparent rate density. The output from this code is processed to
make the figures.

## Extra goodies:

`src/numerical_models.py` implements Monte Carlo surveys for (A) twin binary,
varying planet populations, and for (B) fixed planet radius, varying binary
mass ratio populations. The routine behind this is explained in more depth in
`monte_carlo_readme.md`, and has only been tested for these two cases.
