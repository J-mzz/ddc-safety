function out = ddc_safety_robust(Cons,Region,eps,epsw,options)
% DDC_SAFETY_ROBUST  SOS program for the robust data-driven safe controller.
%   out = ddc_safety_robust(Cons, Region, eps, epsw, options) finds
%   polynomials rho (density) and psi such that the rational controller
%   u = psi/rho keeps every system in the consistency set away from the
%   unsafe set, robustly to process noise ||w||_inf <= epsw (epsw = 0 gives
%   the non-robust controller).
%
%   Cons    problem data: n, T, degrees drho, dpsi, df, dg, and the
%           consistency set N*theta <= e (fields N, e)
%   Region  initial set (r0, c0) and unsafe sets (ru1, cu1), (ru2, cu2);
%           r is the squared radius and c the center of each disk
%   eps     noise bound on the data (the consistency set already uses it)
%
%   out.cr, out.cp  coefficients of rho and psi in the monomial basis
%   out.sol         YALMIP diagnostics (out.F, out.Q, out.v: constraints,
%                   Gram matrices and monomials)

kkk = 1;

A   = Cons.A;
B   = Cons.B;
xi  = Cons.xi;
n   = Cons.n;
T   = Cons.T;
drho = Cons.drho;
dpsi = Cons.dpsi;
df  = Cons.df;
dg  = Cons.dg;
N = Cons.N;
e = Cons.e;

r0 = Region.r0;
c0 = Region.c0;
ru1 = Region.ru1;
cu1 = Region.cu1;
ru2 = Region.ru2;
cu2 = Region.cu2;

% initialization
sdpvar x1 x2
vars = [x1,x2];

[rho,cr] = polynomial(vars,drho,0);
[psi,cp] = polynomial(vars,dpsi,kkk);

[~,~,vf] = polynomial(vars,df,kkk);
[~,~,vg] = polynomial(vars,dg,0);

[s1,c1] = polynomial(vars,drho-2);
[s2,c2] = polynomial(vars,drho-2);

% slakness
[~,~,v] = polynomial(vars, drho/2, 0);
tol1 =  v'*1e-8*eye(length(v))*v;       % slackness, as mu in (20) 

% define initial and unsafe region
x0 = r0 - (x1-c0(1))^2 - (x2-c0(2))^2;	% rho > 0
xu1 = (ru1 - (x1-cu1(1))^2 - (x2-cu1(2))^2);	% rho < 0
xu2 = (ru2 - (x1-cu2(1))^2 - (x2-cu2(2))^2);   % rho < 0
xu_com = -xu1 * xu2;

% solve with yalmip
dmax = max(dpsi+dg-kkk, drho+df-kkk);

Y = zeros(1, size(N,1), 'like', sdpvar);      % the polynomial matrix Y
for i = 1:size(N,1)
    Y(i) = polynomial(vars, dmax, 0);	% see (20) for size
end

YY = zeros(1,2*n,'like',sdpvar);        % extra part for noise
for j = 1:2*n
    YY(j) = polynomial(vars, dmax, 0);
end

% (15) in the note
D1 = [kron(jacobian(rho,vars),vf') kron(jacobian(psi,vars),vg')];
D2 = [rho*[jacobian(vf,x1)' jacobian(vf,x2)'] psi*[jacobian(vg,x1)' jacobian(vg,x2)']];

k = coefficients(-D1-D2-Y*N, vars);

% the slack constants 1e-5 and 10*tol1 are the ones used for Figs. 1 and 2
F = [k == 0,...
    sos(-Y*e - rho*xu_com+1e-5),...
    (sos(Y)):'y',...
    sos(rho -s1*x0 -1*tol1),...	% rho > 0 on X0 (rho = a/b where b > 0)
    sos(-rho-s2*xu_com -10*tol1),...    % rho < 0 on Xu
    sos(-rho*xu_com -psi),...    % |psi| <= -rho*xu_com
    sos(-rho*xu_com +psi),...
    sos(s1),sos(s2),...
    sum(cr) == 1e-7,...
    ];

if epsw ~= 0
    kk = coefficients(-jacobian(rho,vars) -YY*[eye(n);-eye(n)], vars);

    F = [F,...
        kk == 0,...
        sos(-YY*[epsw*ones(n, 1); epsw*ones(n, 1)]),...
        sos(YY)];
end

[sol,v,Q] = solvesos(F, [], options, [coefficients(Y,vars);coefficients(YY,vars);c1;c2;cr;cp]);

out.cr = value(cr);
out.cp = value(cp);
out.rho = rho;
out.sol = sol;
out.F = F;
out.Q = Q;
out.v = v;
out.Y = value(coefficients(Y,vars));
out.D1 = value(coefficients(D1,vars));
out.D2 = value(coefficients(D2,vars));

end