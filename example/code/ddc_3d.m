close all;clear;clc
yalmip('clear')

rng(20)
%% initialization
d_r = 4;


%% system
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
%
%% setup
n 	 = 3;   % # of states
T    = 80;  % # of samples
drho = 4;   % degree of rho
dpsi = 4;   % degree of psi
df   = 3;   % degree of f
dg   = 0;   % degree of g

% noise
eps = 1;    % noise level for robust data driven
epsw = 1;   % noise level of disturbance


%% generate data
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


%% consistency set
%   X: data matrix, columns represent [u,t,x1,x2,x3,dx1,dx2,dx3]

lf = nchoosek(df+n, n) -1;         % length of monomial vector for f(x)
lg = nchoosek(dg+n, n);         % ... for g(x)

syms x1 x2 x3
vars = [x1,x2,x3];

p = monomials(vars, 1:df);      % phi
g = monomials(vars, 0:dg);      % gamma

p_fun = matlabFunction(p);
g_fun = matlabFunction(g);

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

%% ddc safety control design
options = sdpsettings('solver', 'mosek', 'verbose', 0);

sdpvar x1 x2 x3
vars = [x1,x2,x3];

% set x0, xu
R0 = 0.01;  C0 = [-.6; 0; 0.1];        % radius,center
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

% k = coefficients(-D1-D2-Y*[A B; -A -B], vars);
k = coefficients(-D1-D2-Y*N, vars);

ubound = 1; % |psi/rho| <= ubound 

F = [k == 0,...
    sos(-Y*e - rho*xu + 1e-4),... %2e-6; 
    (sos(Y)):'y',...
    sos(rho -s1*x0 -1*tol1),...	% since rho = a/b where b > 0
    sos(-rho-s2*xu -1*tol1),...
%     sos(ubound*rho -psi),...    % bound input u
%     sos(ubound*rho +psi),...
    sos(-rho*xu -psi),...    % bound input u
    sos(-rho*xu +psi),...    
    sos(s1),sos(s2),...
    sum(cr) == 1e-7,...
    ];

if epsw ~= 0
    kk = coefficients(-jacobian(rho,vars) -YY*[eye(n);-eye(n)], vars);

    F = [F,...
        kk == 0,...
        sos(-YY*[epsw*ones(n, 1); epsw*ones(n, 1)]),... % +2*tol1),...
        sos(YY)];
end

sol = solvesos(F, [], options, [coefficients(Y,vars);coefficients(YY,vars);c1;c2;cr;cp])
cr = value(cr);
cp = value(cp);

Y = value(coefficients(Y,vars));
D1 = value(coefficients(D1,vars));
D2 = value(coefficients(D2,vars));


% save('all_data')
% 
% load('all_data')
%% extract solutions
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
vec_psi = monomials(vars, 1:dpsi); % 1:Cons.dpsi
rho = cr'*vec_rho;
psi = cp'*vec_psi;
u = psi/rho;

z0 = R0 - (z1-C0(1))^2 - (z2-C0(2))^2 - (z3-C0(3))^2;   % rho > 0
zu = Ru - (z1-Cu(1))^2 - (z2-Cu(2))^2 - (z3-Cu(3))^2;   % rho < 0

v = monomials(vars,0:drho/2);
tol = v'*1e-8*eye(length(v))*v;       % slackness, as mu in (20) 

rhof = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2)+jacobian(f(3),z3));
psig = jacobian(psi,vars)*g + psi*(jacobian(g(1),z1)+jacobian(g(2),z2)+jacobian(g(3),z3));
div = rhof + psig - rho*zu; %2*tol;

%% plot 3d rho, x0, xu
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

%% plot 3d phase portrait
% len = -.8:0.2:0.5;
% 
% [X,Y,Z] = meshgrid(len,len,len);
% 
% U =   - 2*X.^3 + X/2 - 2*Z.^3 + (5*Z)/2 + Y;
% V = 2*Y.^3 - (3*Y)/2 + 2*Z.^3 - (5*Z)/2 - X;
% W =   2*X.^3 - (3*X)/2 + 2*Y.^3 - Y/2 - 2*Z;
% 
% f4 = quiver3(X,Y,Z,U,V,W,'Color',[0.5,0.5,0.5]);

% hold off

%% plot 3d trajectories
y0 = [-.4;0;0];
F = matlabFunction(f);
U = matlabFunction(u);
Xdot = matlabFunction(f+g*u);

%% 
T = 2;
mu = 0.1;
ode_opt = odeset('RelTol', 1e-6,'AbsTol', 1e-8);
n_traj = 30;

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
    
    while t_all < T
        tmax_curr = exprnd(mu);
        tmax_curr = min(t_all + tmax_curr, T) - t_all;

        noise_curr = eps*(2*rand(3,1)-1);
        
        % open loop
%         if i == 1
            f_curr = @(t, y) F(y(1),y(2),y(3));
%         else 
%             f_curr = @(t, y) F(y(1),y(2),y(3)) + noise_curr;
%         end
        
        % close loop
%         if i == 1
%             f_curr = @(t, y) Xdot(y(1),y(2),y(3));
%         else 
%             f_curr = @(t, y) Xdot(y(1),y(2),y(3)) + noise_curr;
%         end
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
    
%     if i == 1
%         plot3(Ydata(1,:),Ydata(2,:),Ydata(3,:),'r','LineWidth',1)
%     else
        plot3(Ydata(1,:),Ydata(2,:),Ydata(3,:),'b')
%     end
    
    DATA.Traj{i} = Ydata;
    DATA.Input{i} = Udata;
    DATA.Time{i} = Tlog;
    DATA.Switch{i} = switch_times;
end
hold off










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% previous simulation method
%{ 
% timer terminate intergrator after 300s
ode_options6 = odeset('RelTol', 1e-7,'AbsTol', 1e-8, 'MaxStep', 1e-6, 'Events',@mytimer);
ode_options5 = odeset('RelTol', 1e-7,'AbsTol', 1e-8, 'MaxStep', 1e-5, 'Events',@mytimer);
tspan = 0:5;

% [t,y] = ode15s(Ft,tspan,y0,ode_options);
% plot3(y(:,1),y(:,2),y(:,3),'r')
global DATALOG
global data_index
global trajectory_index

DATALOG = [];
trajectory_index = 1;

%% trajecories of system without control with noise
% data_index = 1
% tic
% [t,y] = ode15s(@(t,y) mypoly3(t,y,0,0,F(y(1),y(2),y(3))),tspan,y0,ode_options5);
% % YY1 = DATALOG{trajectory_index}(:, 3);
% % YY2 = DATALOG{trajectory_index}(:, 4);
% % YY3 = DATALOG{trajectory_index}(:, 5);
% 
% for trajectory_index = 2:11
%     data_index = 1
%     tic
%     [t,y] = ode15s(@(t,y) mypoly3(t,y,epsw,0,F(y(1),y(2),y(3)) ),tspan,y0,ode_options5);
% %     YY1 = DATALOG{trajectory_index}(:, 3);
% %     YY2 = DATALOG{trajectory_index}(:, 4);
% %     YY3 = DATALOG{trajectory_index}(:, 5);
% %     plot3(YY1(1:end-1),YY2(1:end-1),YY3(1:end-1),'b') 
% end

%% trajecories of system with control with noise
data_index = 1
tic
[t,y] = ode15s(@(t,y) mypoly3(t,y,0,U(y(1),y(2),y(3)),F(y(1),y(2),y(3))),tspan,y0,ode_options5);

for trajectory_index = 2:31
    data_index = 1
    tic
    [t,y] = ode15s(@(t,y) mypoly3(t,y,epsw,U(y(1),y(2),y(3)),F(y(1),y(2),y(3)) ),tspan,y0,ode_options5);
    YY1 = DATALOG{trajectory_index}(:, 3);
    YY2 = DATALOG{trajectory_index}(:, 4);
    YY3 = DATALOG{trajectory_index}(:, 5);
    plot3(YY1(1:end-1),YY2(1:end-1),YY3(1:end-1),'b') 
end


save('datalog_close_clean1_noise','DATALOG')
%% 

function dxdt = mypoly3(t,x,epsw,u,f)
global DATALOG
global data_index
global trajectory_index

    DATALOG{trajectory_index}(data_index, 1) = u;
    DATALOG{trajectory_index}(data_index, 2) = t;
    DATALOG{trajectory_index}(data_index, 3) = x(1);
    DATALOG{trajectory_index}(data_index, 4) = x(2);
    DATALOG{trajectory_index}(data_index, 5) = x(3);
    noise = epsw*(2*rand(3,1)-1);
    dxdt = f + [0;0;1]*u + noise;
    DATALOG{trajectory_index}(data_index, 6) = dxdt(1);
    DATALOG{trajectory_index}(data_index, 7) = dxdt(2);
    DATALOG{trajectory_index}(data_index, 8) = dxdt(3);
    DATALOG{trajectory_index}(data_index,  9) = noise(1);
    DATALOG{trajectory_index}(data_index, 10) = noise(2);
    DATALOG{trajectory_index}(data_index, 11) = noise(3);
    data_index = data_index + 1;
end



function [value, isterminal, direction] = mytimer(t, y)
% set a timer to terminate in 'TimeOut' seconds
    TimeOut = 300;
    value = toc-TimeOut;
    isterminal = 1;
    direction = 0;
end
% YY1 = DATALOG{1}(:, 3);YY2 = DATALOG{1}(:, 4);YY3 = DATALOG{1}(:, 5);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%