'''
calculate the ratios of apparent to true occurrence rates for the
order-of-magnitude toy model.
'''
nb_by_ns = 0.05

Λa_by_Λ_at_r = (1 + 2**(3/2)*nb_by_ns)**(-1)

Λa_by_Λ_at_rbyrt2 = 2**(5/2)*nb_by_ns / (1 + 2**(3/2)*nb_by_ns)

print(Λa_by_Λ_at_r, Λa_by_Λ_at_r*0.5)
print(Λa_by_Λ_at_rbyrt2, Λa_by_Λ_at_rbyrt2*0.5)

bf = 1/(1+1/nb_by_ns)
print(bf)
