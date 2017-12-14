#!/usr/bin/env bash

# "Nominal" models (this just means I arbitrarily decided Z_2 = 0.5)
python integrate_for_apparent_rate_density.py --modelnumber 1 --ZsubTwo 0.5 --binaryfrac 0.1 --upperradiuscutoff 1
python integrate_for_apparent_rate_density.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
python integrate_for_apparent_rate_density.py --modelnumber 4 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
python integrate_for_apparent_rate_density.py --modelnumber 5 --ZsubTwo 0.5 --binaryfrac 0.1 --upperradiuscutoff 22.5
python integrate_for_apparent_rate_density.py --modelnumber 6 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5

# Check vs Z_2
python integrate_for_apparent_rate_density.py --modelnumber 3 --ZsubTwo 0.0 --binaryfrac 0.44 --upperradiuscutoff 22.5
python integrate_for_apparent_rate_density.py --modelnumber 3 --ZsubTwo 0.25 --binaryfrac 0.44 --upperradiuscutoff 22.5

# Check vs upper radius cutoff
python integrate_for_apparent_rate_density.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 15.0
python integrate_for_apparent_rate_density.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 17.5
python integrate_for_apparent_rate_density.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 20.0
python integrate_for_apparent_rate_density.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 25.0

# If you fine-tune both r_pu AND Z_2/Z_0 preferentially, how big of a "HJ
# discrepancy" do you get?
python integrate_for_apparent_rate_density.py --modelnumber 3 --ZsubTwo 0.0 --binaryfrac 0.44 --upperradiuscutoff 15.0
