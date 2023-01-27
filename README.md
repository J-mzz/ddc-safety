# ddc-safety

Run ddc_nonauto_robust.m, it solves the problem when we consider noise in divergence, ie, div( rho(f+gu+noise) )

epsw is the noise level for that noise, epsw = 0 will reduce the problem to no noise in divergence, ie, div( rho(f+gu) )
eps is the noise level for data-driven method.

(un)comment line 78-79 in ddc_safety_robust.m to (un)bound input
    

For robust case without boundness on the input, for epsw=1e-2, the green line (unsafe) will cross the red contour (rho=0).
                with boundness on the input, this code will return junks.
