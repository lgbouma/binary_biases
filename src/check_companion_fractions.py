from __future__ import division
import numpy as np
from scipy.integrate import trapz

α = 3.5
β = 0
δ = -2.92

q = np.arange(0,1+1e-6,1e-6)

A = (1 + q**α)**(-1/2)

int1 = trapz( q**β * A**(-δ-4), q )
int2 = trapz( q**β * A**(-4), q )

BF = 0.44

frac = (1 + (BF/(1-BF)*int1)**(-1)) / (1 + (BF/(1-BF)*int2)**(-1) )

print(frac)
