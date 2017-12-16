'''
each Monte Carlo iteration asks the Q: how many HJs get detected?

then we ask: what is the probability of getting at least 10 detected HJs if you
survey a sample of 836 stars, given a rate drawn from the stated bounds of
Petigura et al (in prep)?

Made by pythonizing Josh Winn's IDL code asking the same question.
'''
from __future__ import division
import numpy as np

n_stars = 836 # number of stars in Wright et al (2012)
n_HJ = 10 # number of HJs in Wright et al (2012)

m = int(1e5)  # number of monte carlo iterations

eta = 5.7e-3 # Petigura+ in prep's reported HJ rate
eta_unc = 1.3e-3

####################
ndet = np.zeros((m))

seed = 42

for i in range(m):
    np.random.seed(seed+i)
    Lambda = eta + eta_unc*np.random.normal()
    ndet[i] = np.random.poisson(n_stars*Lambda)

q = ndet >= n_HJ

print('p means fraction of MC trials with at least {:d} HJs'.
        format(n_HJ))
print('without taking Z or binaries into account, p = {:.3f}'.format(
    np.sum(q)/m
    ))


####################
eta *= 1.2
ndet = np.zeros((m))

seed = 42

for i in range(m):
    np.random.seed(seed+i)
    Lambda = eta + eta_unc*np.random.normal()
    ndet[i] = np.random.poisson(n_stars*Lambda)

q = ndet >= n_HJ
print('taking Z into account, p = {:.3f}'.format(
    np.sum(q)/m
    ))

####################
eta *= 1.13
ndet = np.zeros((m))

seed = 42

for i in range(m):
    np.random.seed(seed+i)
    Lambda = eta + eta_unc*np.random.normal()
    ndet[i] = np.random.poisson(n_stars*Lambda)

q = ndet >= n_HJ
print('taking Z and binaries into account, p = {:.3f}'.format(
    np.sum(q)/m
    ))
