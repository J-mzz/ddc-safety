function [u] = min_norm_scalar(z,rhof, rhog)
%MIN_NORM_SCALAR
%
%u(z) = argmin_u u^2 | rhof + u rhog >= 0 (tol?)
%
rhof_eval = rhof(z(1), z(2));
rhog_eval = rhog(z(1), z(2));

tol = 1e-6;

if rhof_eval>tol
    u = 0;
else
    u = (tol-rhof_eval)/rhog_eval;
end

div_out = rhof_eval + rhog_eval*u;
end

