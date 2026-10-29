function plot_example2(panel, file)
% PLOT_EXAMPLE2  Draw a panel of Fig. 3 (Example 2, twist system).
%   plot_example2(panel) draws 'fig3a' or 'fig3b' from the results reported
%   in the paper (data/*.mat).
%   plot_example2(panel, file) draws the results saved by example2_twist.m
%   (set save_file there) instead.
%
%   The figure shows the trajectories, X0 (black), Xu (red) and, for the
%   closed loop, the surface rho = 0 (green).

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'utils')))

switch panel
    case 'fig3a'
        data_file = 'fig3a_open_loop.mat';
        name = 'Example 2 open-loop no process noise';
        show_rho = false;
    case 'fig3b'
        data_file = 'fig3b_robust.mat';
        name = 'Example 2 robust with process noise';
        show_rho = true;
    otherwise
        error('Unknown panel ''%s''.', panel)
end
if nargin < 2
    file = fullfile(root, 'data', data_file);
end

S = load(file, 'drho', 'cr', 'R0', 'C0', 'Ru', 'Cu', 'DATA');
R0 = S.R0; C0 = S.C0;
Ru = S.Ru; Cu = S.Cu;

%% rho
syms z1 z2 z3
vars = [z1; z2; z3];

vec_rho = monomials(vars, 0:S.drho);
rho = S.cr'*vec_rho;

%% plot 3d rho, x0, xu
figure()
colormaps = linspecer(6);
hold on
grid on

if show_rho
    f1 = fimplicit3(rho,'FaceColor',colormaps(3,:),'EdgeColor','none');
end
[Xs, Ys, Zs] = sphere(50);
f2 = surf(sqrt(R0)*Xs+C0(1), sqrt(R0)*Ys+C0(2), sqrt(R0)*Zs+C0(3), 'FaceColor','k','EdgeColor','none');
f3 = surf(sqrt(Ru)*Xs+Cu(1), sqrt(Ru)*Ys+Cu(2), sqrt(Ru)*Zs+Cu(3), 'FaceColor','r','EdgeColor','none');
view(3)

xlim([-.8,0.2])
ylim([-.5,0.5])
zlim([-.5,0.5])

%% trajectories
for ind = 1:numel(S.DATA.Traj)
    YY1 = S.DATA.Traj{ind}(1, :);
    YY2 = S.DATA.Traj{ind}(2, :);
    YY3 = S.DATA.Traj{ind}(3, :);
    f5 = plot3(YY1,YY2,YY3,'Color',[0 0.4470 0.7410],'LineWidth',1);
end

hold off

caz = 185.0709;
cel = 67.8552;
view(caz,cel)

if show_rho
    legend([f1,f2,f3,f5],{'\rho = 0','x0','xu','x(t)'},'FontSize',16)
else
    legend([f2,f3,f5],{'x0','xu','x(t)'},'FontSize',16)
end
title(name,'FontSize', 16)

end
