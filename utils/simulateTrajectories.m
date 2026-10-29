function DATA = simulateTrajectories(Region,f,g,u,epsw,sim_mode,n,T)
% SIMULATETRAJECTORIES  Simulate Example 1 from the boundary of X0 and plot.
%   DATA = simulateTrajectories(Region, f, g, u, epsw, sim_mode, n, T)
%   simulates n trajectories on t in [0, T] that start on the boundary of
%   the initial set (center Region.c0, squared radius Region.r0):
%       sim_mode = 'closed':  x_dot = f(x) + g*u(x) + w(t)
%       sim_mode = 'open':    x_dot = f(x) + w(t)
%   The noise w(t) is piecewise constant, uniform in [-epsw, epsw]^2, and
%   switches at exponentially distributed times (mean 0.1). The trajectories
%   and the closed-loop phase portrait are plotted in the current figure.
%   f, g and u are symbolic in z1, z2.
%   Defaults: sim_mode = 'closed', n = 30, T = 2.
%
%   DATA.Traj{i}    states (2 x K)            DATA.Time{i}    time stamps
%   DATA.Input{i}   logged inputs             DATA.Switch{i}  switching times

if nargin < 6, sim_mode = 'closed'; end
if nargin < 7, n = 30; end
if nargin < 8, T = 2; end

r0 = Region.r0;
c0 = Region.c0;     % initial

FF = matlabFunction(f+g*u);
F = matlabFunction(f);
U = matlabFunction(u);

mu = 0.1;
ode_opt = odeset('RelTol', 1e-6,'AbsTol', 1e-8);

% initial states on the boundary of X0
theta = (1:n)/n * 2*pi;
xxx = sqrt(r0)*sin(theta)+c0(1);
yyy = sqrt(r0)*cos(theta)+c0(2);

for i = 1:n
    yprev = [xxx(i);yyy(i)];
    t_all = 0;
    switch_times = 0;

    Ydata = yprev;
    Udata = U(yprev(1),yprev(2));
    Tlog = 0;

    fprintf('current traj: %d \n', i)

    while t_all < T
        tmax_curr = exprnd(mu);
        tmax_curr = min(t_all + tmax_curr, T) - t_all;

        noise_curr = epsw*(2*rand(2,1)-1);

        switch sim_mode
            case 'closed'
                f_curr = @(t, y) FF(y(1),y(2)) + noise_curr;
            case 'open'
                f_curr = @(t, y) F(y(1),y(2)) + noise_curr;
            otherwise
                error('sim_mode must be ''closed'' or ''open''.')
        end

        [tcurr, ycurr] = ode15s(f_curr, [0, tmax_curr], yprev, ode_opt);

        Ydata = [Ydata, ycurr'];
        Udata = [Udata, U(ycurr(1),ycurr(2))];
        Tlog = [Tlog; t_all + tcurr];
        switch_times = [switch_times; t_all + tmax_curr];

        yprev = ycurr(end, :)';
        t_all = t_all + tmax_curr;
    end

    plot(Ydata(1,2:end),Ydata(2,2:end),'b')

    DATA.Traj{i} = Ydata;
    DATA.Input{i} = Udata;
    DATA.Time{i} = Tlog;
    DATA.Switch{i} = switch_times;
end

xlim([-5 5])
ylim([-5 5])
axis square

%% plot phase portrait
plotpp(@(t,x)  F(x(1),x(2)) + g * U(x(1),x(2)))

end
