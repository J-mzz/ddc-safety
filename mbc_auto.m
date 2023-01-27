close all;clear;clc

%% safety
mode = 'barrier'; % 'barrier' or 'density'

yalmip('clear')
% define a,b,c,s1,s2
sdpvar x1 x2
vars = [x1; x2];
d_r = 4;
d_f = 1;

% without integration condition
[s0,c0] = polynomial(vars, d_r-2);
[s1,c1] = polynomial(vars, d_r-2);
[s2,c2] = polynomial(vars, d_r-2);
[s3,c3] = polynomial(vars, d_r-2);

[rho,cr] = polynomial(vars, d_r, 0);      % a

% system
A = [-0.6409    0.8962
    0.2408   -0.6790];
f = A*vars;
f = [x2; -x1 + 1/3*x1^3 - x2];

% set x0, xu
R0 = .25;  C0 = [2; 2];        % radius,center
Ru = .25;  Cu = [-2; -2];
Re = .09;  Ce = [0; 0];

x0 = R0 - (x1-C0(1))^2 - (x2-C0(2))^2;   % rho > 0
xu = Ru - (x1-Cu(1))^2 - (x2-Cu(2))^2;   % rho < 0
xe = Re - (x1-Ce(1))^2 - (x2-Ce(2))^2;

% compute div rho
if strcmp(mode,'barrier')
    divrho = jacobian(rho,vars)*f;
elseif strcmp(mode,'density')
    divrho = jacobian(rho,vars)*f + rho*(jacobian(f(1),x1)+jacobian(f(2),x2));
end

% constraints
F = [sos(divrho + xe*s0),... %  - .6
     sos(rho - x0*s1), ...
     sos(-rho - xu*s2), ...
     sos(rho - xe*s3), ...
     sos(s1),sos(s2),sos(s0),sos(s3)];
options = sdpsettings('solver','mosek','verbose',0);
sol = solvesos(F, [], options, [cr;c1;c2;c0;c3])

% extract solutions for plot
cr = value(cr);
syms z1 z2
vars = [z1;z2];
vec_rho = monomials(vars, 0:d_r);
rho = cr'*vec_rho;

% system
f = A*vars;
f = [z2; -z1 + 1/3*z1^3 - z2];

if strcmp(mode,'barrier')
    divrho = jacobian(rho,vars)*f;
elseif strcmp(mode,'density')
    divrho = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2));
end

% plot
f1 = fcontour(rho,'r');
f1.LevelList = 0;
hold on
f2 = sdisplay(x0);
f2 = f2{1};
f2 = str2sym(f2);
f2 = fcontour(f2,'k');
f2.LevelList = 0;
f3 = sdisplay(xu);
f3 = f3{1};
f3 = str2sym(f3);
f3 = fcontour(f3,'g');
f3.LevelList = 0;

f4 = fcontour(divrho,'m');
f4.LevelList= 0;

f5 = sdisplay(xe);
f5 = f5{1};
f5 = str2sym(f5);
f5 = fcontour(f5,'b');
f5.LevelList = 0;

RHO = matlabFunction(rho);

%{

% compute rho, drho, divrho
divrho = jacobian(rho,vars)*f + rho*jacobian(f,x1)+jacobian(f,x2);

RHO = matlabFunction(rho);
DIVRHO = matlabFunction(divrho);

% phase portrait
q = 5;     % bound
k = 0.2;   % step size
[X,Y] = meshgrid(-q:k:q, -q:k:q);

for i = 1:size(X,1)
    for j = 1:size(X,2)
    U(i,j) = A(1,:)*[X(i,j);Y(i,j)];
    V(i,j) = A(2,:)*[X(i,j);Y(i,j)];
    end
end
% V = V + subs(u,{x1 x2}, {X Y});

figure
quiver(X,Y,U,V, 'autoscalefactor', 10)
[startx, starty] = meshgrid(-q:q:q,-q:q:q);
streamline(X,Y,U,V,startx,starty, 0.2)
axis(1*[-q q -q q])
xlabel('x_1', 'fontsize', 13)
ylabel('x_2', 'fontsize', 13)
set(gca,'fontsize', 13)
legend('Phase Portrait','FontSize',15)

%}



% plot_traj_safety
% test_plot


