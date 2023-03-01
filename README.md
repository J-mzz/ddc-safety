# ddc-safety

ddc_nonauto_robust.m solves the robust problem, ie, div( rho(f+gu+noise) )

epsw is the online noise level, epsw = 0 will reduce to nonrobust case.
eps is the offline noise level for data-driven method.

Current slackness setting for figures.

Robsut:     sos (- Ye + 3.00E-05) & sos (-rho-s2*xu - 8.8*tol1)

Nonrobust:  sos (- Ye + 1.00E-05) & sos (-rho-s2*xu - 10 *tol1)

Open loop:  sos (- Ye + 1.00E-05) & sos (-rho-s2*xu - 10 *tol1)

