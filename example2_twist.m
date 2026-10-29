%% Example 2 (twist system) of "Robust Data-Driven Safe Control Using Density Functions"
%
% True system (only used to generate data and to simulate):
%   x_dot = Af*x + Bf*(4*x.^3 - 3*x)/2 + [0; 0; 1]*u + w,   ||w||_inf <= epsw
% Initial set X0: ball centered at (-0.6, 0, 0.1) with radius 0.1
% Unsafe set Xu:  ball centered at (-0.1, 0, 0.1) with radius 0.1
% Prior: f is a cubic polynomial with f(0) = 0 and g is constant; 80 samples
%        of x_dot in [-0.8, 0.8]^3 with noise ||.||_inf <= eps = 1
%
% The script solves the SOS program for the density rho and psi (control
% u = psi/rho), plots rho = 0, X0 and Xu, and simulates 30 trajectories
% that start on a circle on the boundary of X0.
%
% Set `panel` below to the panel of Fig. 3 to reproduce. The figures in the
% paper are drawn with plot_example2.m from data/.

close all;clear;clc
yalmip('clear')

%% Settings
% panel of the paper:
%   'fig3a'  open loop, no process noise
%   'fig3b'  robust controller, with process noise ||w||_inf <= 1
panel = 'fig3b';

eps  = 1;   % noise bound on the data of x_dot
epsw = 1;   % process noise bound in the controller design

% sim_mode: 'open' (u = 0) or 'closed' (u = psi/rho) in the simulation
% epsw_sim: process noise bound in the simulation
% T_sim:    simulation horizon
switch panel
    case 'fig3a', sim_mode = 'open';   epsw_sim = 0; T_sim = 5;
    case 'fig3b', sim_mode = 'closed'; epsw_sim = 1; T_sim = 5;
    otherwise
        error('Unknown panel ''%s''.', panel)
end

% file to save the results for plot_example2(panel, save_file); '' = do not save
save_file = '';

rng(20)
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')))

%% System
sdpvar x1 x2 x3
vars = [x1; x2; x3];

Af = [-1  1  1;
     -1  0 -1;
      0  1 -2];
Bf = [-1  0 -1;
      0  1  1;
      1  1  0];

temp_vars = 1/2*[4*x1^3 - 3*x1;
                 4*x2^3 - 3*x2;
                 4*x3^3 - 3*x3];

f = Af*vars + Bf*temp_vars;
g = [0; 0; 1];

%% Problem dimensions
n 	 = 3;   % # of states
T    = 80;  % # of samples
drho = 4;   % degree of rho
dpsi = 4;   % degree of psi
df   = 3;   % degree of f
dg   = 0;   % degree of g

%% Generate data
box = .8;

syms z1 z2 z3
vars = [z1;z2;z3];

temp_vars = 1/2*[4*z1^3 - 3*z1;
                 4*z2^3 - 3*z2;
                 4*z3^3 - 3*z3];

f = Af*vars + Bf*temp_vars;

sample_data = [];
for i = 1:T

    u = 2*rand(1,1)-1;

	noise = eps*(2*rand(n,1)-1);

    xsample = (2*rand(n,1)-1) * box;

    z1 = xsample(1);
    z2 = xsample(2);
    z3 = xsample(3);

    xdot = subs(f+g*u+noise);

    curr = [u,i,z1,z2,z3,xdot(1),xdot(2),xdot(3)];

    sample_data = [sample_data;curr];
end

sample_data = double(sample_data);


%% Consistency set
%   sample_data: one sample per row [u,t,x1,x2,x3,dx1,dx2,dx3]

lf = nchoosek(df+n, n) -1;      % length of monomial vector for f(x)
lg = nchoosek(dg+n, n);         % ... for g(x)

syms x1 x2 x3
vars = [x1,x2,x3];

p = monomials(vars, 1:df);      % phi
gam = monomials(vars, 0:dg);    % gamma

p_fun = matlabFunction(p);
g_fun = matlabFunction(gam);

A = []; B = [];
p_data  = zeros(T, lf);
ug_data = zeros(T, lg);

ind = randperm(size(sample_data,1),T);
sample_data = sample_data(ind,:);

for i = 1:T
    p_data(i,:) = feval(p_fun, sample_data(i,3), sample_data(i,4), sample_data(i,5));

    if dg == 0
        ug_data(i,:) = sample_data(i, 1);
    elseif dg > 0
        ug_data(i, :) = sample_data(i, 1)*feval(g_fun, sample_data(i,3), sample_data(i,4), sample_data(i,5));
    end

    A = [A; kron(eye(n),p_data(i,:))];
    B = [B; kron(eye(n),ug_data(i,:))];
end

xi = sample_data(1:T,[6 7 8])';
xi = xi(:);

% nontrivial constraints
N = [A B; -A -B];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end

e = [eps*ones(n*T,1)+xi; eps*ones(n*T,1)-xi];
[N_red, e_red] = nontrivial_constraints(N, e);
N = N_red;
e = e_red;

%% Data-driven safe control design (SOS program)
options = sdpsettings('solver', 'mosek', 'verbose', 0);

sdpvar x1 x2 x3
vars = [x1,x2,x3];

% set x0, xu
R0 = 0.01;  C0 = [-.6; 0; 0.1];        % squared radius, center
Ru = 0.01;  Cu = [-.1; 0; 0.1];

x0 = R0 - (x1-C0(1))^2 - (x2-C0(2))^2 - (x3-C0(3))^2;   % rho > 0
xu = Ru - (x1-Cu(1))^2 - (x2-Cu(2))^2 - (x3-Cu(3))^2;   % rho < 0

% initialization
[rho,cr] = polynomial(vars,drho,0);
[psi,cp] = polynomial(vars,dpsi,1);

[~,~,vf] = polynomial(vars,df,1);
[~,~,vg] = polynomial(vars,dg,0);

[s1,c1] = polynomial(vars,drho-2);
[s2,c2] = polynomial(vars,drho-2);

% slackness
[~,~,v] = polynomial(vars, drho/2, 0);
tol1 =  v'*1e-8*eye(length(v))*v;       % slackness, as mu in (20)

% solve with yalmip
dmax = max(dpsi+dg-1, drho+df-1);

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
D2 = [rho*[jacobian(vf,x1)' jacobian(vf,x2)' jacobian(vf,x3)'] ...
    psi*[jacobian(vg,x1)' jacobian(vg,x2)' jacobian(vg,x3)']];

k = coefficients(-D1-D2-Y*N, vars);

F = [k == 0,...
    sos(-Y*e - rho*xu + 1e-4),...
    (sos(Y)):'y',...
    sos(rho -s1*x0 -1*tol1),...	% since rho = a/b where b > 0
    sos(-rho-s2*xu -1*tol1),...
    sos(-rho*xu -psi),...    % bound input u
    sos(-rho*xu +psi),...
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

sol = solvesos(F, [], options, [coefficients(Y,vars);coefficients(YY,vars);c1;c2;cr;cp])
cr = value(cr);
cp = value(cp);

Y = value(coefficients(Y,vars));
D1 = value(coefficients(D1,vars));
D2 = value(coefficients(D2,vars));

%% Extract the solution: rho, psi, and u = psi/rho
syms z1 z2 z3
vars = [z1; z2; z3];

Af = [-1  1  1;
     -1  0 -1;
      0  1 -2];
Bf = [-1  0 -1;
      0  1  1;
      1  1  0];

temp_vars = 1/2*[4*z1^3 - 3*z1;
                 4*z2^3 - 3*z2;
                 4*z3^3 - 3*z3];

f = Af*vars + Bf*temp_vars;
g = [0; 0; 1];

vec_rho = monomials(vars, 0:drho);
vec_psi = monomials(vars, 1:dpsi);
rho = cr'*vec_rho;
psi = cp'*vec_psi;
u = psi/rho;

z0 = R0 - (z1-C0(1))^2 - (z2-C0(2))^2 - (z3-C0(3))^2;   % rho > 0
zu = Ru - (z1-Cu(1))^2 - (z2-Cu(2))^2 - (z3-Cu(3))^2;   % rho < 0

v = monomials(vars,0:drho/2);
tol = v'*1e-8*eye(length(v))*v;       % slackness, as mu in (20)

rhof = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2)+jacobian(f(3),z3));
psig = jacobian(psi,vars)*g + psi*(jacobian(g(1),z1)+jacobian(g(2),z2)+jacobian(g(3),z3));
div = rhof + psig - rho*zu;

%% Plot rho = 0, X0 and Xu
figure()
f1 = fimplicit3(rho,'FaceColor',[0.3010 0.7450 0.9330],'EdgeColor','none');
hold on

[Xs, Ys, Zs] = sphere(50);
surf(sqrt(R0)*Xs+C0(1), sqrt(R0)*Ys+C0(2), sqrt(R0)*Zs+C0(3), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none')
surf(sqrt(Ru)*Xs+Cu(1), sqrt(Ru)*Ys+Cu(2), sqrt(Ru)*Zs+Cu(3), 'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none')
view(3)

xlim([-.85,0.15])
ylim([-.5,0.5])
zlim([-.5,0.5])

%% Simulate trajectories
F = matlabFunction(f);
U = matlabFunction(u);
Xdot = matlabFunction(f+g*u);

mu = 0.1;   % mean time between switches of the noise
ode_opt = odeset('RelTol', 1e-6,'AbsTol', 1e-8);
n_traj = 30;

% initial states on a circle on the boundary of X0
theta = (1:n_traj)/n_traj * 2*pi;
xxx = -0.6+.1*cos(theta);
yyy = 0.1*sin(theta);
zzz = 0.1;

for i = 1:n_traj
    yprev = [xxx(i);yyy(i);zzz];
    t_all = 0;
    switch_times = 0;

    Ydata = yprev;
    Udata = U(yprev(1),yprev(2),yprev(3));
    Tlog = 0;

    fprintf('current traj: %d \n', i)

    while t_all < T_sim
        tmax_curr = exprnd(mu);
        tmax_curr = min(t_all + tmax_curr, T_sim) - t_all;

        noise_curr = epsw_sim*(2*rand(3,1)-1);

        switch sim_mode
            case 'open'
                f_curr = @(t, y) F(y(1),y(2),y(3)) + noise_curr;
            case 'closed'
                f_curr = @(t, y) Xdot(y(1),y(2),y(3)) + noise_curr;
        end

        warning off
        [tcurr, ycurr] = ode15s(f_curr, [0, tmax_curr], yprev, ode_opt);
        warning on

        Ydata = [Ydata, ycurr'];
        Udata = [Udata, U(ycurr(1),ycurr(2),ycurr(3))];
        Tlog = [Tlog; t_all + tcurr];
        switch_times = [switch_times; t_all + tmax_curr];

        yprev = ycurr(end, :)';
        t_all = t_all + tmax_curr;
    end

    plot3(Ydata(1,:),Ydata(2,:),Ydata(3,:),'b')

    DATA.Traj{i} = Ydata;
    DATA.Input{i} = Udata;
    DATA.Time{i} = Tlog;
    DATA.Switch{i} = switch_times;
end
hold off

%% Save the results
if ~isempty(save_file)
    save(save_file, 'panel', 'eps', 'epsw', 'epsw_sim', 'sim_mode', 'T_sim', ...
        'n', 'drho', 'dpsi', 'R0', 'C0', 'Ru', 'Cu', 'cr', 'cp', 'DATA')
end
