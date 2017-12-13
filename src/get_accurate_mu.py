import numpy as np
from scipy.integrate import trapz

α = 3.5
β = 0

q = np.arange(0,1+1e-8,1e-8)

# 2**(3/2) + 3I = this
this = trapz( (1 + q**α)**(3/2), q )

I = (this - 2**(3/2))/3

print(-I)
