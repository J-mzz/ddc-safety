function x = getDataRND(eg,eps,N)
% GETDATARND  N noisy samples of x_dot at random states.
%   x = getDataRND(eg, eps, N) samples x in [-5, 5]^2 and u in [-1, 1] and
%   returns one row [u, i, x1, x2, dx1, dx2] per sample, where
%   dx = f(x) + g*u + w and ||w||_inf <= eps (see getSystem.m).

box = 5;

syms x1 x2
vars = [x1;x2];
[f,g] = getSystem(eg,vars);

x = [];
for i = 1:N
    
    u = 2*rand(1,1)-1;

	noise = eps*(2*rand(2,1)-1);
    
    x0 = (2*rand(2,1)-1) * box;
    
    x1 = x0(1);
    x2 = x0(2);
    
    xdot = subs(f+g*u+noise);
    
    curr = [u,i,x1,x2,xdot(1),xdot(2)];
    
    x = [x;curr];
end

x = double(x);
end