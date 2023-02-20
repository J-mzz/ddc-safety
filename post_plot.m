close all;clear;clc
load 15_w_u_same_x0_7

%x0 = [DATA.X{1,i}(1); DATA.Y{1,i}(1)];

R0 = Region.r0;
C0 = Region.c0;     % initial
Ru = Region.ru;
Cu = Region.cu;
x0 = R0 - (x1-C0(1))^2 - (x2-C0(2))^2;   % rho > 0
xu = Ru - (x1-Cu(1))^2 - (x2-Cu(2))^2;   % rho < 0

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
div = rhof + psig + 2e-6;

R0 = Region.r0;
C0 = Region.c0;     % initial
Ru = Region.ru;
Cu = Region.cu;
x0 = R0 - (x1-C0(1))^2 - (x2-C0(2))^2;   % rho > 0
xu = Ru - (x1-Cu(1))^2 - (x2-Cu(2))^2;   % rho < 0

%% plot
figure(1)
clf
f1 = fcontour(rho,'r','LineWidth',1);
f1.LevelList = 0;
hold on
f2 = sdisplay(x0);
f2 = f2{1};
f2 = str2sym(f2);
f2 = fcontour(f2,'k','LineWidth',1);
f2.LevelList = 0;
f3 = sdisplay(xu);
f3 = f3{1};
f3 = str2sym(f3);
f3 = fcontour(f3,'g','LineWidth',1);
f3.LevelList = 0;
% f4 = fcontour(div,'b');
% f4.LevelList = 0;

% check function value
DIV = matlabFunction(div);
RHO = matlabFunction(rho);
PSI = matlabFunction(psi);
U = matlabFunction(psi/rho);
FF = matlabFunction(f+g*u);
F = matlabFunction(f);

%% plot phase portrait and trajectories for closed loop

warning off

n = 10; % # of samples
sample = (1:n)/n * 2*pi;
x = C0(1) + sqrt(R0) * sin(sample);
y = C0(2) + sqrt(R0) * cos(sample);

for i = 1:n

    xdata = DATA.X{i};
    ydata = DATA.Y{i};
    f5 = plot(xdata, ydata,'Color',[0.3010 0.7450 0.9330],'LineWidth',1);
    
end

xxx0 = DATA.X{1,1}(1)
yyy0 = DATA.Y{1,1}(1)

plot(xxx0, yyy0,'bo','LineWidth',1)
% text(xxx0, yyy0,['(' num2str(xxx0) ',' num2str(yyy0) ')'])

%% plot phase portrait

plotpp(@(t,x)  F(x(1),x(2)) + g * U(x(1),x(2)))

% legend([f1,f2,f3,f4],{'rho','x0','xu','div'},'FontSize',12)
legend([f1,f2,f3,f5],{'rho','x0','xu','x(t)'},'FontSize',12)
warning on


xlim([-5 5])
ylim([-5 5])
axis square