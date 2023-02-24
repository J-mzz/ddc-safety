close all;clear;clc
% load data_2d_open_2unsafe
load temp2


sdpvar x1 x2
syms z1 z2
vars = [z1; z2];
[f,g] = getSystem(eg,vars);	% symbolic

vec_rho = monomials(vars, 0:Cons.drho);
vec_psi = monomials(vars, 1:Cons.dpsi); % 1:Cons.dpsi
rho = cr'*vec_rho;
psi = cp'*vec_psi;
u = psi/rho;

v = monomials(vars,0:Cons.drho/2);
tol = v'*1e-8*eye(length(v))*v;       % slackness, as mu in (20) 

rhof = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2));
psig = jacobian(psi,vars)*g + psi*(jacobian(g(1),z1)+jacobian(g(2),z2));
% div = rhof + psig + 2e-6;

R0 = Region.r0;
C0 = Region.c0;     % initial
Ru = Region.ru;
Cu = Region.cu;
x0 = R0 - (x1-C0(1))^2 - (x2-C0(2))^2;   % rho > 0
xu = Ru - (x1-Cu(1))^2 - (x2-Cu(2))^2;   % rho < 0
xu1 = 0.16 - (x1+1)^2 - (x2-1)^2;   % rho < 0

%% plot
figure(1)
clf
hold on
colormaps = linspecer(6);
%% plot trajectories

warning off

n = length(DATA.Traj);

for i = 1:n
%     if i == 1
%         data = DATA.Traj{i};
%         f5 = plot(data(1,:), data(2,:), 'Color',[0 0.4470 0.7410],'LineWidth',2); %colormaps(5,:)
%     else
        data = DATA.Traj{i};
        f5 = plot(data(1,:), data(2,:), 'Color',[0.3010 0.7450 0.9330],'LineWidth',.8);% colormaps(2,:)
%     end
end

data0 = data(:,1);

% plot(data0(1), data0(2),'bo','LineWidth',1)

%% plot level set

f1 = fcontour(rho,'LineColor',colormaps(3,:),'LineWidth',1); %'g'
f1.LevelList = 0;
f2 = sdisplay(x0);
f2 = f2{1};
f2 = str2sym(f2);
f2 = fcontour(f2,'k','LineWidth',1); % 'LineColor',colormaps(4,:)
f2.LevelList = 0;
f3 = sdisplay(xu);
f3 = f3{1};
f3 = str2sym(f3);
f3 = fcontour(f3,'r','LineWidth',1); % 'LineColor',colormaps(1,:)
f3.LevelList = 0;
% f4 = fcontour(div,'b');
% f4.LevelList = 0;
f6 = sdisplay(xu1);
f6 = f6{1};
f6 = str2sym(f6);
f6 = fcontour(f6,'r','LineWidth',1); % 'LineColor',colormaps(1,:)
f6.LevelList = 0;

% check function value
DIV = matlabFunction(div);
RHO = matlabFunction(rho);
PSI = matlabFunction(psi);
U = matlabFunction(psi/rho);
FF = matlabFunction(f+g*u);
F = matlabFunction(f);

%% plot phase portrait

plotpp(@(t,x)  F(x(1),x(2)) + g * U(x(1),x(2)))

% legend([f1,f2,f3,f4],{'rho','x0','xu','div'},'FontSize',12)
legend([f1,f2,f3,f5],{'rho','x0','xu','x(t)'},'FontSize',12)
warning on
title('Example 1 closed-loop with 2 unsafe regions')

xlim([-5.5 3.5])
ylim([-5 4])
axis square