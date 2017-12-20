import pandas as pd, numpy as np
import matplotlib.pyplot as plt

#from NASA exoplanet archive, only KOIs with "CONFIRMED" designation
df = pd.read_csv('../data/confirmed_KOI_radii.csv', comment='#')

plt.hist(np.log10(df['koi_prad']), bins=20)
plt.savefig('../results/KOI_radii_log.pdf', bbox_inches='tight')

print(len(df[df['koi_prad']>2])/len(df))
