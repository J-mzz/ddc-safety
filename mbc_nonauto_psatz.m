close all;clear;clc

%% safety
mode = 'density'; % 'barrier' or 'density'

yalmip('clear')
% define a,b,c,s1,s2
sdpvar x1 x2
vars = [x1; x2];
d = 6;
% d=4;

x = [x1; x2];
% system

g = [0; 1];

% f = [x2; -x1 - x2];

f = [x2; -x1 + 1/3*x1^3 - x2];

[rho,cr] = polynomial(vars, d, 0);      % rho
[psi,cp] = polynomial(vars, d, 0);      % psi = rho*u

% slackness
% [~,~,v] = polynomial(vars, d/2, 0);  
% tol1 =  v'*1e-8*eye(length(v))*v;       % slackness in rho

tol1 = 1e-5;
tolset = 1e-3;

% set x0, xu
R0 = 0.25;  C0 = [0; -3];        % radius,center
Ru = 0.16;  Cu = [-1; -1];
Re = 0.09;  Ce = [0; 0];

x0 = R0 - (x1-C0(1))^2 - (x2-C0(2))^2;   % rho > 0
xu = Ru - (x1-Cu(1))^2 - (x2-Cu(2))^2;   % rho < 0
xe = Re - (x1-Ce(1))^2 - (x2-Ce(2))^2;

% compute div rho
if strcmp(mode,'barrier')
    div = jacobian(rho,vars)*f + ...
             jacobian(psi,vars)*g;
elseif strcmp(mode,'density')
    div = jacobian(rho,vars)*f + rho*(jacobian(f(1),x1)+jacobian(f(2),x2)) + ...
             jacobian(psi,vars)*g + psi*(jacobian(g(1),x1)+jacobian(g(2),x2));
end

X0 = struct('ineq', x0, 'eq', []);
Xu = struct('ineq', [], 'eq', xu);
Xall = struct('ineq', -xu, 'eq', []);
[p_div, con_div, coeff_div] = constraint_psatz(div-tol1, Xall, x, d+2);
[p_init, con_init, coeff_init] = constraint_psatz(rho -tolset, X0, x, d);
[p_unsafe, con_unsafe, coeff_unsafe] = constraint_psatz(-rho -tolset, Xu, x, d);


BOUND_PSI = false;

if BOUND_PSI
    [p_upos, con_upos, coeff_upos] = constraint_psatz(-xu*rho + psi, Xall, x, d);
    [p_uneg, con_uneg, coeff_uneg] = constraint_psatz(-xu*rho - psi, Xall, x, d);
else
    con_upos = [];
    con_uneg = [];
    coeff_upos = [];
    coeff_uneg = [];
end

F = [con_div; con_init; con_unsafe; con_upos; con_uneg];
coeff = [cr; cp; coeff_div; coeff_init; coeff_unsafe; coeff_upos;coeff_uneg];

% constraints
% F = [sos(div - tol1), ...
%      sos(rho  - x0*s1 -1e-3), ...     % r >= 0, in initial x0
%      sos(-rho - xu*s2 -1e-3), ...     % r <= 0, in unsafe xu
%      sos(-xu*rho + psi), ...
%      sos(-xu*rho - psi), ...
%      sos(s1)];

% [p_out, cons, coeff_list] = constraint_psatz(p, X, vars, d)
% sos(s2) in constraint, only need to be <0 on boundary of x0 if starting
% in a different region
options = sdpsettings('solver','mosek','verbose',2);
sol = solvesos(F, norm(cp, 'inf'), options, coeff)

%%  extract solutions for plot
cr = value(cr);
cp = value(cp);
nanind = find(isnan(cp));
cp(nanind) = 0;
cp;

syms z1 z2
vars = [z1;z2];
vec_rho = monomials(vars, 0:d);
vec_psi = monomials(vars, 0:d);
rho = cr'*vec_rho;
psi = cp'*vec_psi;
u = psi/rho;
U = matlabFunction(u);

% v = monomials(vars,0:d/2);
% tol =  v'*1e-8*eye(length(v))*v;       % slackness, as mu in (20) 
% tol = 1e-5;
tol = tol1;
% system
f = [z2; -z1 + 1/3*z1^3 - z2];

if strcmp(mode,'barrier')
    div = jacobian(rho,vars)*f + ...
       jacobian(psi,vars)*g;
elseif strcmp(mode,'density')
    div = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2)) + ...
        jacobian(psi,vars)*g + psi*(jacobian(g(1),z1)+jacobian(g(2),z2))+2*tol;
end

%% plot
% contour line and (un)safe regiobs
f1 = fcontour(rho,'r');
f1.LevelList = 0;
hold on
f2 = sdisplay(x0);
f2 = f2{1};
f2 = str2sym(f2);
f2 = fcontour(f2,'k');
f2.LevelList = 0;
f3 = sdisplay(xu);
f3 = f3{1};
f3 = str2sym(f3);
f3 = fcontour(f3,'g');
f3.LevelList = 0;
% f4 = fcontour(div,'b');
% f4.LevelList = 0;

DIV = matlabFunction(div);
RHO = matlabFunction(rho);
PSI = matlabFunction(psi);

Region.r0 = 0.25;  Region.c0 = [0; -3];        % radius,center
Region.cu = 0.16;  Region.cu = [-1; -1];
% Re = 0.09;  Ce = [0; 0];

% phase portrait and closed loop dynamics
warning off
DATA_rational = testSys(Region,f,g,u,1e-3);
warning on


