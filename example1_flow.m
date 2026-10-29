%% Example 1 (flow system) of "Robust Data-Driven Safe Control Using Density Functions"
%
% True system (only used to generate data and to simulate):
%   x_dot = [x2; -x1 + 1/3*x1^3 - x2] + [0; 1]*u + w,   ||w||_inf <= epsw
% Initial set X0: disk centered at (0, -3) with radius 0.5
% Unsafe set Xu:  disks centered at (-1, -1) and (-1, 1) with radius 0.4
% Prior: f is a cubic polynomial with f(0) = 0 and g is constant; 80 samples
%        of x_dot with noise ||.||_inf <= eps = 2
%
% The script solves the SOS program for the density rho and psi (control
% u = psi/rho), plots the level sets rho = 0, X0, Xu, and simulates
% trajectories from the boundary of X0 with the closed-loop phase portrait.
%
% Set `panel` below to the panel of Figs. 1 and 2 to reproduce. The
% figures in the paper are drawn with plot_example1.m from data/.

close all;clear;clc
yalmip('clear')

%% Settings
% panel of the paper:
%   'fig1a'  open loop, no process noise
%   'fig1b'  robust controller (designed for ||w||_inf <= 2), with process noise
%   'fig2a'  non-robust controller (designed for w = 0), no process noise
%   'fig2b'  non-robust controller, with process noise ||w||_inf <= 2
panel = 'fig1b';

eps = 2;    % noise bound on the data of x_dot

% epsw:     process noise bound in the controller design (0 = non-robust)
% sim_mode: 'open' (u = 0) or 'closed' (u = psi/rho) in the simulation
% epsw_sim: process noise bound in the simulation
% n_traj:   number of simulated trajectories
% T_sim:    simulation horizon
switch panel
    case 'fig1a', epsw = 0; sim_mode = 'open';   epsw_sim = 0; n_traj = 30;  T_sim = 5;
    case 'fig1b', epsw = 2; sim_mode = 'closed'; epsw_sim = 2; n_traj = 30;  T_sim = 2;
    case 'fig2a', epsw = 0; sim_mode = 'closed'; epsw_sim = 0; n_traj = 30;  T_sim = 2;
    case 'fig2b', epsw = 0; sim_mode = 'closed'; epsw_sim = 2; n_traj = 100; T_sim = 2;
    otherwise
        error('Unknown panel ''%s''.', panel)
end

% file to save the results for plot_example1(panel, save_file); '' = do not save
save_file = '';

rng(10)
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')))

%% Initial and unsafe sets (squared radius, center)
Region.r0 = 0.25;
Region.c0 = [0; -3];     % initial

Region.ru1 = 0.16;
Region.cu1 = [-1; -1];   % unsafe

Region.ru2 = 0.16;
Region.cu2 = [-1; 1];    % unsafe

%% System and data
% x_dot = f(x) + g(x)*u(x) + noise, see utils/getSystem.m
eg = 2;

sdpvar x1 x2
vars_sdp = [x1; x2];
[f,g] = getSystem(eg,vars_sdp);   % sdpvar

% problem dimensions and degrees
Cons = setCons(eg);

% data, with input u \in [-1,1]
X = getDataRND(eg,eps,Cons.T);

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

%% Data-driven safe control design (SOS program)
options = sdpsettings('solver', 'mosek', 'verbose', 1);
out = ddc_safety_robust(Cons,Region,eps,epsw,options);

out.sol
cr = out.cr;
cp = out.cp;

%% Extract the solution: rho, psi, and u = psi/rho
syms z1 z2
vars = [z1; z2];
[f,g] = getSystem(eg,vars);	% symbolic

vec_rho = monomials(vars, 0:Cons.drho);
vec_psi = monomials(vars, 1:Cons.dpsi);
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

% divergence of rho*(f + g*u)
rhof = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2));
psig = jacobian(psi,vars)*g + psi*(jacobian(g(1),z1)+jacobian(g(2),z2));
div = rhof + psig;

%% Plot rho = 0, X0 and Xu
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

% function handles, e.g. to check function values
DIV = matlabFunction(div);
RHO = matlabFunction(rho);
PSI = matlabFunction(psi);
U = matlabFunction(u);
DotX = matlabFunction(f+g*u);

xlim([-5,5])
ylim([-5,5])
axis square
grid on

%% Simulate trajectories and plot the closed-loop phase portrait
warning off
DATA = simulateTrajectories(Region,f,g,u,epsw_sim,sim_mode,n_traj,T_sim);

legend([f1,f2,f3],{'rho','x0','xu'},'FontSize',12)
warning on

%% Save the results
if ~isempty(save_file)
    save(save_file, 'panel', 'eg', 'eps', 'epsw', 'epsw_sim', 'sim_mode', ...
        'T_sim', 'Cons', 'Region', 'cr', 'cp', 'DATA')
end
