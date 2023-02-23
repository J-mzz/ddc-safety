close all;clear;clc
yalmip('clear')

%% safety

% define a,b,c,s1,s2
sdpvar x1 x2 x3
vars = [x1; x2; x3];

d_r = 4;

% system
A = [-1  1  1;
     -1  0 -1;
      0  1 -2];
B = [-1  0 -1;
      0  1  1;
      1  1  0];

temp_vars = 1/2*[4*x1^3 - 3*x1;
                 4*x2^3 - 3*x2;
                 4*x3^3 - 3*x3];

f = A*vars + B*temp_vars;
g = [0; 0; 1];

% set x0, xu
R0 = 0.04;  C0 = [-.5; 0; 0];        % radius,center
Ru = 0.04;  Cu = [0; 0; 0];

x0 = R0 - (x1-C0(1))^2 - (x2-C0(2))^2 - (x3-C0(3))^2;   % rho > 0
xu = Ru - (x1-Cu(1))^2 - (x2-Cu(2))^2 - (x3-Cu(3))^2;   % rho < 0

% without integration condition
[s1,c1] = polynomial(vars, d_r-2);
[s2,c2] = polynomial(vars, d_r-2);

[rho,cr] = polynomial(vars, d_r, 0);      % rho
[psi,cp] = polynomial(vars, d_r, 0);      % psi = rho*u

% slackness
[~,~,v] = polynomial(vars, d_r/2, 0);  
tol1 =  v'*1e-8*eye(length(v))*v;       % slackness in rho

% compute div rho
div = jacobian(rho,vars)*f + ...
    rho*(jacobian(f(1),x1) + jacobian(f(2),x2) + jacobian(f(3),x3)) + ...
    jacobian(psi,vars)*g + ...
    psi*(jacobian(g(1),x1) + jacobian(g(2),x2) + jacobian(g(3),x3));

ubound = 1; % bound on input: |psi/rho| < ubound

% constraints
F = [sos(div + 1e-6), ...
     sos(rho  - x0*s1), ...     % r >= 0, in initial x0
     sos(-rho - xu*s2 -1e-8), ...     % r <= 0, in unsafe xu
     sos(ubound*rho + psi), ...
     sos(ubound*rho - psi), ...
     sos(s1),sos(s2)];
options = sdpsettings('solver','mosek','verbose',0);
sol = solvesos(F, [], options, [cr;c1;c2;cp])

% extract solutions for plot
cr = value(cr);
cp = value(cp);
nanind = find(isnan(cp));
cp(nanind) = 0;
cp;

syms z1 z2 z3
vars = [z1; z2; z3];
vec_rho = monomials(vars, 0:d_r);
vec_psi = monomials(vars, 0:d_r);
rho = cr'*vec_rho;
psi = cp'*vec_psi;
u = psi/rho;

% v = monomials(vars,0:d_r/2);
% tol =  v'*1e-8*eye(length(v))*v;

% system
temp_vars = 1/2*[4*z1^3 - 3*z1;
                 4*z2^3 - 3*z2;
                 4*z3^3 - 3*z3];

f = A*vars + B*temp_vars;

z0 = R0 - (z1-C0(1))^2 - (z2-C0(2))^2 - (z3-C0(3))^2;   % rho > 0
zu = Ru - (z1-Cu(1))^2 - (z2-Cu(2))^2 - (z3-Cu(3))^2;   % rho < 0

div = jacobian(rho,vars)*f + ...
    rho*(jacobian(f(1),z1) + jacobian(f(2),z2) + jacobian(f(3),z3)) + ...
    jacobian(psi,vars)*g + ...
    psi*(jacobian(g(1),z1) + jacobian(g(2),z2) + jacobian(g(3),z3));

DIV = matlabFunction(div);
RHO = matlabFunction(rho);
PSI = matlabFunction(psi);
U = matlabFunction(u);



%% plot 3d rho, x0, xu
figure()
f1 = fimplicit3(rho,'FaceColor',[0.3010 0.7450 0.9330],'EdgeColor','none');
hold on 
f2 = fimplicit3(z0,'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none');
f3 = fimplicit3(zu,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none');

xlim([-.8,0.5])
ylim([-.8,0.5])
zlim([-.8,0.5])

%% plot 3d phase portrait
len = -.8:0.2:0.5;

[X,Y,Z] = meshgrid(len,len,len);

U =   - 2*X.^3 + X/2 - 2*Z.^3 + (5*Z)/2 + Y;
V = 2*Y.^3 - (3*Y)/2 + 2*Z.^3 - (5*Z)/2 - X;
W =   2*X.^3 - (3*X)/2 + 2*Y.^3 - Y/2 - 2*Z;

f4 = quiver3(X,Y,Z,U,V,W,'Color',[0.5,0.5,0.5]);


%% plot 3d trajectories
y0 = [-.3;0;0];
F = matlabFunction(f);
Ft = @(t,y) F(y(1),y(2),y(3));

ode_options = odeset('RelTol', 1e-7,'AbsTol', 1e-8, 'MaxStep', 1e-6);
tspan = 0:5;
[t,y] = ode45(Ft,tspan,y0,ode_options);

plot3(y(:,1),y(:,2),y(:,3),'r')
hold off