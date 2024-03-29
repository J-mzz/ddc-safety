close all;clear;clc
yalmip('clear')

rng(10)
%% safe/ unsafe region
% radius, center
Region.r0 = 0.25;
% Region.c0 = [1.5; 0];     % initial
Region.c0 = [0; -3];     % initial
% Region.c0 = [1; 0];     % initial

Region.ru1 = 0.16;
Region.cu1 = [-1; -1];     % unsafe

Region.ru2 = 0.16;
Region.cu2 = [-1; 1];
%% generate system and data
% define system: x_dot = f(x) + g(x)*u(x) + noise

eg = 2; % 1 (linear) or 2 (nonlinear)

sdpvar x1 x2
vars_sdp = [x1; x2];
[f,g] = getSystem(eg,vars_sdp);   % sdpvar

% setup
Cons = setCons(eg);

% data, with input u \in [-1,1]
eps = 2; % noise level for robust data driven
epsw = 2; % noise level of disturbance
%% even for unbounded control: eps 1e-3 works but 1e-2 not 

X = getDataRND(eg,eps,Cons.T); %

% consistency set
[A,B,xi] = getCons(X,Cons);

N = [A B; -A -B];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end

e = [eps*ones(Cons.n*Cons.T,1)+xi; eps*ones(Cons.n*Cons.T,1)-xi];
[N_red, e_red] = nontrivial_constraints(N, e);
Cons.N = N_red;
Cons.e = e_red;

Cons.A = A;
Cons.B = B;
Cons.xi = xi;

%% ddc safety control design
options = sdpsettings('solver', 'mosek', 'verbose', 1);
out = ddc_safety_robust(Cons,Region,eps,epsw,options);

out.sol
cr = out.cr;
cp = out.cp;

% for jjj = 1:size(out.Q,2)
% eigval(jjj) = min(eig(out.Q{jjj}));
% end
% min(eigval)

%% extract solutions
syms z1 z2
vars = [z1; z2];
[f,g] = getSystem(eg,vars);	% symbolic

vec_rho = monomials(vars, 0:Cons.drho);
vec_psi = monomials(vars, 1:Cons.dpsi); % 1:Cons.dpsi
rho = cr'*vec_rho;
psi = cp'*vec_psi;
u = psi/rho;

R0 = Region.r0;
C0 = Region.c0;     % initial
Ru1 = Region.ru1;
Cu1 = Region.cu1;
Ru2 = Region.ru2;
Cu2 = Region.cu2;

x0 = R0 - (x1-C0(1))^2 - (x2-C0(2))^2;   % rho > 0
xu1 = (Ru1 - (x1-Cu1(1))^2 - (x2-Cu1(2))^2);   % rho < 0
xu2 = (Ru2 - (x1-Cu2(1))^2 - (x2-Cu2(2))^2);   % rho < 0
xu_com = -(Ru1 - (x1-Cu1(1))^2 - (x2-Cu1(2))^2) * (Ru2 - (x1-Cu2(1))^2 - (x2-Cu2(2))^2);

zu_com = -(Ru1 - (z1-Cu1(1))^2 - (z2-Cu1(2))^2) * (Ru2 - (z1-Cu2(1))^2 - (z2-Cu2(2))^2);

rhof = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2));
psig = jacobian(psi,vars)*g + psi*(jacobian(g(1),z1)+jacobian(g(2),z2));
div = rhof + psig;

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
f3 = sdisplay(xu_com);
f3 = f3{1};
f3 = str2sym(f3);
f3 = fcontour(f3,'g');
f3.LevelList = 0;
% f4 = fcontour(div,'b');
% f4.LevelList = 0;
% f5 = fcontour(rho*zu,'c');
% f6 = sdisplay(xu1);
% f6 = f6{1};
% f6 = str2sym(f6);
% f6 = fcontour(f6,'g');
% f6.LevelList = 0;

% check function value
DIV = matlabFunction(div);
RHO = matlabFunction(rho);
PSI = matlabFunction(psi);
U = matlabFunction(u);
DotX = matlabFunction(f+g*u);

xlim([-5,5])
ylim([-5,5])
axis square
grid on
%% plot phase portrait and trajectories for closed loop

warning off
DATA = testSys(Region,f,g,u,epsw);

% DATA = New2_testSys(Region,f,g,u,epsw,rho,psi);

legend([f1,f2,f3],{'rho','x0','xu'},'FontSize',12)
warning on

% figure(2)
% clf
% fsurf(rho)
% view(3)

% out.sol.info

% save('9_w_u_same_x0')

% data_X = DATA.X{1};
% data_Y = DATA.Y{1};
% data_U = U(data_X,data_Y);
% data_DX = diff(data_X);
% data_DY = diff(data_Y);
% data_noise = [data_DX;data_DY] - DotX(data_X(1:end-1),data_Y(1:end-1));
% 
% [minimum,index] = min(DATA.U);
% good = [1 5 9 15 19 30]
% bad = [15 16 22]

% for jjj = 1:30
%     jjj
%     min(DATA.Traj{jjj},[],2)
% end