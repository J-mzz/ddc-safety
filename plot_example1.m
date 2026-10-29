function plot_example1(panel, file)
% PLOT_EXAMPLE1  Draw a panel of Figs. 1 and 2 (Example 1, flow system).
%   plot_example1(panel) draws 'fig1a', 'fig1b', 'fig2a' or 'fig2b' from the
%   results reported in the paper (data/*.mat).
%   plot_example1(panel, file) draws the results saved by example1_flow.m
%   (set save_file there) instead.
%
%   The figure shows the closed-loop phase portrait, the trajectories, and
%   the sets rho = 0 (green), X0 (black) and Xu (red).

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'utils')))

% panel settings: data file, title, axis limits, trajectories to plot
traj_idx = [];                  % [] = all trajectories
show_rho = true;
show_legend = true;
switch panel
    case 'fig1a'
        data_file = 'fig1a_open_loop.mat';
        name = 'Example 1 open-loop no process noise';
        lims = [-5.5 3.5 -5 4];
        show_rho = false;
    case 'fig1b'
        data_file = 'fig1b_robust.mat';
        name = 'Example 1 robust with process noise';
        lims = [-5.5 3.5 -5 4];
    case 'fig2a'
        data_file = 'fig2a_nonrobust_no_noise.mat';
        name = 'Example 1 nonrobust without process noise';
        lims = [-4.5 2.5 -4 3];
    case 'fig2b'
        data_file = 'fig2b_nonrobust_noise.mat';
        name = 'Example 1 nonrobust with process noise';
        lims = [-2 0 -2 0];     % zoom on the lower unsafe disk
        traj_idx = [1:4:100 3 12 15 31 100];
        show_legend = false;
    otherwise
        error('Unknown panel ''%s''.', panel)
end
if nargin < 2
    file = fullfile(root, 'data', data_file);
end

S = load(file, 'eg', 'Cons', 'cr', 'cp', 'Region', 'DATA');
Cons = S.Cons;
Region = S.Region;
DATA = S.DATA;

%% rho, psi and u = psi/rho
syms z1 z2
vars = [z1; z2];
[f,g] = getSystem(S.eg,vars);	% symbolic

vec_rho = monomials(vars, 0:Cons.drho);
vec_psi = monomials(vars, 1:Cons.dpsi);
rho = S.cr'*vec_rho;
psi = S.cp'*vec_psi;

R0 = Region.r0;
C0 = Region.c0;     % initial
Ru1 = Region.ru1;
Cu1 = Region.cu1;
Ru2 = Region.ru2;
Cu2 = Region.cu2;

x0 = R0 - (z1-C0(1))^2 - (z2-C0(2))^2;   % rho > 0
xu = -(Ru1 - (z1-Cu1(1))^2 - (z2-Cu1(2))^2) * (Ru2 - (z1-Cu2(1))^2 - (z2-Cu2(2))^2);   % rho < 0

U = matlabFunction(psi/rho);
F = matlabFunction(f);

%% plot
figure()
hold on
colormaps = linspecer(6);

% phase portrait of the closed loop
plotpp(@(t,x)  F(x(1),x(2)) + g * U(x(1),x(2)))

% trajectories
if isempty(traj_idx)
    traj_idx = 1:length(DATA.Traj);
end
for i = traj_idx
    data = DATA.Traj{i};
    f5 = plot(data(1,:), data(2,:), 'Color',[0 0.4470 0.7410],'LineWidth',1);
end

% level sets
f1 = fcontour(rho,'LineColor',colormaps(3,:),'LineWidth',2);
f1.LevelList = 0;
f1.Visible = show_rho;
f2 = fcontour(x0,'k','LineWidth',2);
f2.LevelList = 0;
f3 = fcontour(xu,'r','LineWidth',2);
f3.LevelList = 0;

if show_legend && show_rho
    legend([f1,f2,f3,f5],{'\rho = 0','x0','xu','x(t)'},'FontSize',16)
elseif show_legend
    legend([f2,f3,f5],{'x0','xu','x(t)'},'FontSize',16)
end
title(name,'FontSize', 16)

xlim(lims(1:2))
ylim(lims(3:4))
axis square
grid on
hold off

end
