BF = 0.05
nb_by_ns = 1 / ( (1/BF) - 1 )

Λa_by_Λ_at_r = (1 + 2**(3/2)*nb_by_ns)**(-1)

Λa_by_Λ_at_rbyrt2 = 2**(5/2)*nb_by_ns / (1 + 2**(3/2)*nb_by_ns)

print(Λa_by_Λ_at_r, Λa_by_Λ_at_r*0.5)
print(Λa_by_Λ_at_rbyrt2, Λa_by_Λ_at_rbyrt2*0.5)
