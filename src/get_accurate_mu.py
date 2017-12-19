from __future__ import division
import numpy as np
from scipy.integrate import trapz

α = 3.5
β = 0

#q = np.arange(0,1+1e-8,1e-8)
q = np.arange(0,1+1e-7,1e-7)

# 2**(3/2) + 3I = this
this = trapz( q**β * (1 + q**α)**(3/2), q )

I = (this - 2**(3/2))/3

print(-I)

BF = 0.1
mu = BF/(1-BF) * this
print(mu)

δ = -2.92
print( 2**((δ+3)/2) )
print( (1 + 2**((δ+3)/2) * mu)/(1 + mu))
