----------
DESCRIPTION

This program simulates the following toy surveys:

MODEL #1: Fixed stars, fixed planets, twin binaries.

MODEL #2: Fixed planets and primaries, varying secondaries.

MODEL #3: Fixed primaries, varying planets and secondaries.

In the "nominal" models, all classes of errors are included when calculating
what the observer infers. These errors are

1. misinterpretation of planetary radii (b/c of transit depth dillution; also
   b/c of wrongly assumed host star radii),
2. incorrectly assumed number of selected stars.

The detection efficiency, rather surprisingly, does not seem to be biased in
any leading-order manner.

----------
The simulation works as follows.

First, the user specifies their inputs: the binary fraction, the model class
(#1,2,3), various exponents (α,β,γ,δ), and the `Z_i`'s.

The population of selected stars is constructed as follows.  First, we note
that each selected star has a `star_type` (single, primary, secondary), a
binary mass ratio, if applicable, and the property of whether it is
"searchable".

Every star in the main `numerical_transit_survey` function is assumed to be
searchable (for planets of whatever size are assigned). This doesn't make any
sense if the stars are though to be part of some volume-limited sample. This is
explicitly the point -- the only stars being counted are the ones that are
searchable. At a fixed planet radius and period, this corresponds to a
magnitude-limited sample of stars.

The absolute number of stars is arbitrary; we take 10^6 single stars
throughout. The number of binaries is calculated according to analytic
formulae. The binary mass ratios are, when applicable, samples from the
appropriate magnitude-limited distribution (given α and β).

The procedure for assigning planets is then as follows:
    * each selected star of type i gets a planet at rate `Λ_i`
    * the radii of planets are assigned independently of any host system
    property, as sampled from p_r(r) ~ r^δ for Model #3 or a δ function for
    Models #1 and #2.
    * a planet is detected when a) it transits, and and b) its host star is
    searchable.

The probability of transiting single stars in our model is assumed to be known,
and so it can be mostly corrected by the observer attempting to infer an
occurrence rate. The only bias is for secondaries, which can be smaller. This
effect is including when computing the transit probability.

For detected planets, apparent radii are computed according to analytic
formulae that account for both dilution and the misclassification of stellar
radii. We assume that the observer thinks they observe the primary.

The rates are then computed in bins of true planet radius and apparent planet
radius.

In a given radius bin, the true rate is found by counting the number of planets
that exist around selected stars of all types (singles, primaries,
secondaries), and dividing by the total number of these stars.

The apparent rates are found by counting the number of detected planets that
were found in an apparent radius bin, dividing by the geometric transit
probability for single stars, and dividing by the apparent total number of
stars.

----------
USAGE

Change parameters in the "inputs" section below. Run with Python3.X

----------
The following tests/sanity checks are implemented:

* γ == 0
* there are as many selected primaries as secondaries
* stars in the same system get the same mass ratio
* all mass ratios are > 0
* only things with detected planets have an apparent radius
* the sum of the true distribution's histogrammed rate densities is the true
  average occurrence rate, to three decimals

Model #1:
* numerical value of `X_Γ` at `r_p` matches analytic prediction

Model #2:
* numerical value of `X_Γ` at `r_p` matches analytic prediction
* the inferred rate from 0.73rtrue to 0.99rtrue is that predicted analytically,
  to within 0.5%. (No more precise b/c # of Poisson noise in bins).

Model #3:
* the numerically realized TRUE rate (for the selected stars) agrees with the
  analytic one (from 3rearth to 22rearth).

Model #4:
* none

Model #5:
* the analytic and numeric (true, for selected stars) Λ(r) distributions, above
  3rearth and below 22rearth, agree.
* the analytic and numeric (apparent) Λa(ra) distributions, above ra ~=
  2r\oplus and below rpu/sqrt(2) agree to within 1%.
