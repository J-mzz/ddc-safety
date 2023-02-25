close all;clear;clc

load('all_data')
rng(10)
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

xlim([-.8,0.2])
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
T = 1;
mu = 0.01;
ode_opt = odeset('RelTol', 1e-6,'AbsTol', 1e-8);


for i = 1:11
    yprev = y0;
    t_all = 0;
    switch_times = 0;

    Ydata = y0;
    Udata = U(y0(1),y0(2),y0(3));
    Tlog = 0;
    
    fprintf('current traj: %d \n', i)
    
    while t_all < T
        tmax_curr = exprnd(mu);
        tmax_curr = min(t_all + tmax_curr, T) - t_all;

        noise_curr = epsw*(2*rand(3,1)-1);
        
        % open loop
        if i == 1
            f_curr = @(t, y) F(y(1),y(2),y(3));
        else 
            f_curr = @(t, y) F(y(1),y(2),y(3)) + noise_curr;
        end
        
        % close loop
%         if i == 1
%             f_curr = @(t, y) Xdot(y(1),y(2),y(3));
%         else 
%             f_curr = @(t, y) Xdot(y(1),y(2),y(3)) + noise_curr;
%         end
        
        [tcurr, ycurr] = ode45(f_curr, [0, tmax_curr], yprev, ode_opt);

        Ydata = [Ydata, ycurr'];
        Udata = [Udata, U(ycurr(1),ycurr(2),ycurr(3))];
        Tlog = [Tlog; t_all + tcurr];
        switch_times = [switch_times; t_all + tmax_curr];

        yprev = ycurr(end, :)';
        t_all = t_all + tmax_curr;
    end
    
    if i == 1
        plot3(Ydata(1,:),Ydata(2,:),Ydata(3,:),'r','LineWidth',1)
    else
        plot3(Ydata(1,:),Ydata(2,:),Ydata(3,:),'b')
    end
    
    DATA.Traj{i} = Ydata;
    DATA.Input{i} = Udata;
    DATA.Time{i} = Tlog;
    DATA.Switch{i} = switch_times;
end
hold off