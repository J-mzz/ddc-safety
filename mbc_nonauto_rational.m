close all;clear;clc

%% safety
mode = 'density'; % 'barrier' or 'density'

yalmip('clear')
% define a,b,c,s1,s2
sdpvar x1 x2
vars = [x1; x2];
d_r = 4;

% system
A = [-0.6409    0.8962
    0.2408   -0.6790];
g = [0; 1];

f = A*vars;
f = [x2; -x1 + 1/3*x1^3 - x2];
% f = [x2; (1-x1^2)*x2-x1];

% without integration condition
[s0,c0] = polynomial(vars, d_r-2);
[s1,c1] = polynomial(vars, d_r-2);
[s2,c2] = polynomial(vars, d_r-2);
[s3,c3] = polynomial(vars, d_r-2);

[rho,cr] = polynomial(vars, d_r, 0);      % a, rho = a/b
[psi,cp] = polynomial(vars, d_r, 0);      % c, psi = rho*u = c/b

% with integration condition: rho = a/b, psi = c/b
% option 1: fixed b

[~,~,v] = polynomial(vars, 3, 0);
b = v'*eye(length(v))*v;

% option 2: polynomial b, then moment

% [b,cb] = polynomial(vars,d_r,0);



% slackness
[~,~,v] = polynomial(vars, d_r/2, 0);  
tol1 =  v'*1e-8*eye(length(v))*v;       % slackness in rho

% set x0, xu
R0 = 0.25;  C0 = [1.5; 0];        % radius,center
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
    div = (jacobian(rho,vars)*b - jacobian(b,vars)*rho) * f + ...
          rho*b * (jacobian(f(1),x1)+jacobian(f(2),x2)) + ...
          (jacobian(psi,vars)*b - jacobian(b,vars)*psi) * g + ...
          psi*b * (jacobian(g(1),x1)+jacobian(g(2),x2));
end

ubound = 1; % bound on input: |psi/rho| < ubound

% constraints
% without eq point

F = [sos(div + 2*tol1), ...
     sos(rho - x0*s1 -tol1), ...      % r >= 0, in initial x0
     sos(-rho - xu*s2 -tol1), ...     % r <= 0, in unsafe xu
     sos(ubound*rho + psi), ...
     sos(ubound*rho - psi), ...
     sos(s1),sos(s2)];
options = sdpsettings('solver','mosek','verbose',1);
sol = solvesos(F, [], options, [cr;c1;c2;cp])

%% extract solutions for plot
cr = value(cr)
cp = value(cp)
nanind = find(isnan(cp));
cp(nanind) = 0;
cp;

syms z1 z2
vars = [z1;z2];
vec_rho = monomials(vars, 0:d_r);
vec_psi = monomials(vars, 0:d_r);

v = monomials(vars,0:1);
b = v'*eye(length(v))*v;

rho = cr'*vec_rho;
psi = cp'*vec_psi;
u = psi/rho;
U = matlabFunction(u);

v = monomials(vars,0:d_r/2);
tol =  v'*1e-8*eye(length(v))*v;       % slackness, as mu in (20) 

% system
f = A*vars;
f = [z2; -z1 + 1/3*z1^3 - z2];
% f = [z2; (1-z1^2)*z2-z1];

if strcmp(mode,'barrier')
    div = jacobian(rho,vars)*f + ...
       jacobian(psi,vars)*g;
elseif strcmp(mode,'density')
    div = (jacobian(rho,vars)*b - jacobian(b,vars)*rho) * f + ...
          rho*b * (jacobian(f(1),z1)+jacobian(f(2),z2)) + ...
          (jacobian(psi,vars)*b - jacobian(b,vars)*psi) * g + ...
          psi*b * (jacobian(g(1),z1)+jacobian(g(2),z2)) + 2*tol;
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
f4 = fcontour(div,'b');
f4.LevelList = 0;

DIV = matlabFunction(div);
RHO = matlabFunction(rho);

Region.r0 = 0.25;  Region.c0 = [1.5; 0];        % radius,center
Region.cu = 0.16;  Region.cu = [-2; -2];
% Re = 0.09;  Ce = [0; 0];

% phase portrait and closed loop dynamics
warning off
testSys(Region,f,g,u,0)
warning on


