'''
This is a test I wrote that fails because my analytics are wrong.

It would be ideal if we could compare the analytic and numeric (apparent)
Λ_a(r_a) distributions, above r_a ~= 2r_\oplus and below r_pu/sqrt(2).  By
analogy with the fixed planet radius, variable mass ratio model -- this might
seemingly work.  However if you think about it more closely, it won't. You can
derive an analytic result for the f(r)~r^δ case.  But numerically, we have
finite edge effects. A r=20rearth planet can be diluted to much below
20rearth/sqrt(2), if it orbits a secondary.  So to make a good analytic
prediction, you need to include the finite edge effects.  This might be
possible, but it'd certainly be painful.  We've validated enough other cases
that the numerics should be trusted (provided they're tested for a variety of
upper cutoffs on the radius distribution).
'''

del vals_num, vals_anal
Λ_a_num = np.array(inferred_dict['Λ'])

q_grid = np.arange(1e-6, 1+1e-6, 1e-6)
prob_ml_q = q_grid**β * (1+q_grid**α)**(3/2)
norm_q = trapz(prob_ml_q, q_grid)
# Derived 2017/11/30.3 (Eqs 70-72 in email to JNW and KM 2017/11/30).
I_1 = Λ_1/(norm_q*norm_r)*r_a_anal**δ * trapz(
        q_grid**β * (1+q_grid**α)**(δ/2 + 2), q_grid
        )
I_2 = Λ_2/(norm_q*norm_r)*r_a_anal**δ * trapz(
        q_grid**(β+δ+8/3) * (1+q_grid**α)**(3/2) * (1+q_grid**(-α))**((δ+1)/2),
        q_grid
        )
Γ_0 = Λ_0 * r_a_anal**δ / norm_r

Λ_a_anal = 1/(1+μ) * (
        Γ_0 + μ * (I_1+I_2)
        )

vals_num = Λ_a_num[(r_a_anal>3) & (r_a_anal<=r_pu/2**(1/2))]

assert np.allclose(np.diff(inferred_dict['r']),
                   np.diff(inferred_dict['r'])[0])
bin_width = np.diff(true_dict['r'])[0]

vals_anal = bin_width * Λ_a_anal[(r_a_anal>3) & (r_a_anal<=r_pu/2**(1/2))]

if debugging:
    plt.close('all')
    f,ax = plt.subplots()
    ax.plot(r_anal[(r_anal>3) & (r_anal<=r_pu/2**(1/2))], vals_num, label='numeric')
    ax.plot(r_anal[(r_anal>3) & (r_anal<=r_pu/2**(1/2))], vals_anal, label='analytic')
    ax.set_yscale('log')
    ax.legend(loc='best')
    f.savefig('temp_apparent.pdf')

# THE FOLLOWING TEST FAILS: it tries to ensure that the numerically realized
# values for the apparent distribution are within 1% of the expected analytic
# values.
try:
    assert np.isclose(np.mean(vals_num/vals_anal), 1, rtol=1e-2)
except:
    print('hi')
    import IPython; IPython.embed()
print('u made it')
import IPython; IPython.embed()
