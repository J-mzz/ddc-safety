close all;clear;clc

load('robust_3d')

%% extract solutions
syms z1 z2 z3
vars = [z1; z2; z3];

A = [-1  1  1;
     -1  0 -1;
      0  1 -2];
B = [-1  0 -1;
      0  1  1;
      1  1  0];

temp_vars = 1/2*[4*z1^3 - 3*z1;
                 4*z2^3 - 3*z2;
                 4*z3^3 - 3*z3];

f = A*vars + B*temp_vars;
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
colormaps = linspecer(6);
hold on
grid on

f1 = fimplicit3(rho,'FaceColor',colormaps(3,:),'EdgeColor','none');%[0.3010 0.7450 0.9330]
[Xs, Ys, Zs] = sphere(50);
f2 = surf(sqrt(R0)*Xs+C0(1), sqrt(R0)*Ys+C0(2), sqrt(R0)*Zs+C0(3), 'FaceColor','k','EdgeColor','none');%[0.9290 0.6940 0.1250]
f3 = surf(sqrt(Ru)*Xs+Cu(1), sqrt(Ru)*Ys+Cu(2), sqrt(Ru)*Zs+Cu(3), 'FaceColor','r','EdgeColor','none');%[0.4660 0.6740 0.1880]
view(3)

% f2 = fimplicit3(z0,'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none');
% f3 = fimplicit3(zu,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none');

xlim([-.8,0.2])
ylim([-.5,0.5])
zlim([-.5,0.5])

%}
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
Ft = @(t,y) F(y(1),y(2),y(3));
Noise = @(t, y) epsw*(2*rand(n,1)-ones(n,1));

% Ft_noise_open = @(t, y) Xdot_open(y(1),y(2),y(3)) + Noise(t,[y(1),y(2),y(3)]);
Ft_noise = @(t, y) F(y(1),y(2),y(3)) + Noise(t,[y(1),y(2),y(3)]);



% plot open-loop trajectory
% ind = 1;
% YY1 = DATALOG{ind}(:, 3);
% YY2 = DATALOG{ind}(:, 4);
% YY3 = DATALOG{ind}(:, 5);
% plot3(YY1,YY2,YY3,'r')


% load('datalog_step6')

% plot closed-loop trajectory with online noise
for ind = 1:30
    YY1 = DATA.Traj{ind}(1, :);
    YY2 = DATA.Traj{ind}(2, :);
    YY3 = DATA.Traj{ind}(3, :);
    f5 = plot3(YY1,YY2,YY3,'Color',[0 0.4470 0.7410],'LineWidth',1);%[0.3010 0.7450 0.9330]
end

hold off

% caz = 7.8000;
% cel = 14.3772;

% caz = 169.7950;
% cel = 56.4386;

caz = 185.0709;
cel = 67.8552;

view(caz,cel)




legend([f1,f2,f3,f5],{'rho','x0','xu','x(t)'},'FontSize',16)
warning on
title('Example 2 robust with process noise','FontSize', 20)