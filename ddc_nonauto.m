close all;clear;clc
yalmip('clear')

rng(1)
%% safe/ unsafe region
% radius, center
Region.r0 = 0.25;
Region.c0 = [1.5; 0];     % initial

Region.ru = 0.16;
Region.cu = [-1; -1];     % unsafe

%% generate system and data
% define system: x_dot = f(x) + g(x)*u(x) + noise

eg = 2; % 1 (linear) or 2 (nonlinear)

sdpvar x1 x2
vars_sdp = [x1; x2];
[f,g] = getSystem(eg,vars_sdp);   % sdpvar

% setup
Cons = setCons(eg);

% data, with input u \in [-1,1]
eps = 1e-2; % noise level for robust data driven
%% only for bounded control: eps 1e-2 works but 1e-1 not 

X = getDataRND(eg,eps,Cons.T); %

% consistency set
[A,B,xi] = getCons(X,Cons);

% [A_red, b_red] = nontrivial_constraints(A, B);
Cons.A = A;
Cons.B = B;
Cons.xi = xi;

%% ddc safety control design
options = sdpsettings('solver', 'mosek', 'verbose', 0);
out = ddc_safety(Cons,Region,eps,options);

out.sol
cr = out.cr;
cp = out.cp;

%% extract solutions
syms z1 z2
vars = [z1; z2];
[f,g] = getSystem(eg,vars);	% symbolic

vec_rho = monomials(vars, 0:Cons.drho);
vec_psi = monomials(vars, 1:Cons.dpsi);
rho = cr'*vec_rho;
psi = cp'*vec_psi;
u = psi/rho;

rhof = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2));
psig = jacobian(psi,vars)*g + psi*(jacobian(g(1),z1)+jacobian(g(2),z2));
div = rhof + psig;

R0 = Region.r0;
C0 = Region.c0;     % initial
Ru = Region.ru;
Cu = Region.cu;
x0 = R0 - (x1-C0(1))^2 - (x2-C0(2))^2;   % rho > 0
xu = Ru - (x1-Cu(1))^2 - (x2-Cu(2))^2;   % rho < 0

%% plot
figure(1)
clf
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

vec_s = monomials(vars, 0:(Cons.drho-2));
s = out.c1'*vec_s;
vec_t = monomials(vars, 0:(Cons.drho/2));
tol = vec_t'*1e-8*eye(length(vec_t))*vec_t;
x0 = R0 - (z1-C0(1))^2 - (z2-C0(2))^2;   % rho > 0
xu = Ru - (z1-Cu(1))^2 - (z2-Cu(2))^2;   % rho < 0

% check function value
DIV = matlabFunction(div);
RHO = matlabFunction(rho);
PSI = matlabFunction(psi);
U = matlabFunction(psi/rho);
S = matlabFunction(rho-s*x0-tol);

%% plot phase portrait and trajectories for closed loop

warning off
testSys(Region,f,g,u)
warning on

legend([f1,f2,f3,f4],{'rho','x0','xu','div'},'FontSize',12)

% figure(2)
% clf
% fsurf(rho)
% view(3)

% out.sol.info


