%sample from the twist system
%
%uses sampler code from https://github.com/Jarmill/peak/tree/master/utils

%% set up the sampler
s_opt = sampler_options;
s_opt.mu = 0.2;
% s_opt = C0;
BOX = 3;
INIT_POINT = 1;
C0 = [-0.5; 0; 0];        % radius,center
if INIT_POINT
    s_opt.sample.x = @() C0;
end
s_opt.Nb = 3;
s_opt.Tmax = 5;
s_opt.parallel = 0;
dynamics =struct('Tmax', 5, 'discrete', 0, 'time_indep', 0);
dynamics.event = @(t, x) box_event(t, x, BOX);

% Box = 5;
% dynamics.Xval = constraint_func(obj.opts.X, obj.vars.x);
% dynamics.event = {@(t,x,w) support_event(t, x, dynamics.Xval, ...
%       0, obj.opts.Tmax)};


f0 = @(t, x) Xdot(x(1), x(2), x(3)); %dynamics
dynamics.f = @(t,x,w,d,b) f0(t, x) + (2*b-1)*epsw;

%% sample the trajectories
Nsample_traj = 100;
out_sim = sampler(dynamics, Nsample_traj, s_opt);

%% plot
F = figure(30);
clf
hold on 
for i = 1:length(out_sim)
    tcurr = out_sim{i}.t;
    xcurr = out_sim{i}.x;
    plot3(xcurr(:, 1), xcurr(:, 2), xcurr(:, 3), 'c')
end

FS_axis = 14;

xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', FS_axis);
ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', FS_axis);
zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', FS_axis);
% title('Phase Plane', 'FontSize', obj.FS_title);   

% figure()
f1 = fimplicit3(rho,'FaceColor',[0.3010 0.7450 0.9330],'EdgeColor','none');
% hold on 
[Xs, Ys, Zs] = sphere(50);

surf(sqrt(R0)*Xs+C0(1), sqrt(R0)*Ys+C0(2), sqrt(R0)*Zs+C0(3), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none')
surf(sqrt(Ru)*Xs+Cu(1), sqrt(Ru)*Ys+Cu(2), sqrt(Ru)*Zs+Cu(3), 'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none')
view(3)
% f2 = fimplicit3(z0,'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none');
% f3 = fimplicit3(zu,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none');

xlim([-.8,0.5])
ylim([-.8,0.5])
zlim([-.8,0.5])
