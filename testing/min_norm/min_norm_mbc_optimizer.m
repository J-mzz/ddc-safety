function [P] = min_norm_mbc_optimizer(th0, th1, tol)
%get a yalmip optimizer for the data-driven (robust) min-norm controller 
u_min = sdpvar(1);
z = sdpvar(2, 1);

rhof = th0(z(1), z(2));
rhog = th1(z(1), z(2));
% tol1 = tol(z(1), z(2));
tol1 = tol;
% tol1 = 0;

F = [
%     rhof+rhog*u_min + 1e-8 >= 0,...
    rhof+rhog*u_min - tol1 >= 0,...
    ];

options = sdpsettings('solver', 'mosek', 'verbose', 0);
% optimize(F, norm(u_min,2), options)


P = optimizer(F, norm(u_min, 2), options,{z}, u_min);

% u = value(u_min);

end