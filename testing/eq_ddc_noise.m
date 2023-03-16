n = 2;
x = sdpvar(n, 1);
order = 6;
vars = x;
d = 2*order;

TWO_DISK = 1;
SEPARATE_YE = 0;

rng(1)

%set geometry
C0 = [1.5; 0];
R0 = 0.5;

Cu = [0; -1];
Ru = 0.6;

Cu2 = [1; 1];
Ru2 = 0.5;


%support sets
X0 = struct('ineq', R0^2 - sum((x - C0).^2), 'eq', []);
Xu1 = struct('ineq', Ru^2 - sum((x - Cu).^2), 'eq', []);
Xu2 = struct('ineq', Ru2^2 - sum((x - Cu2).^2), 'eq', []);
X = struct('ineq', [sum((x - Cu).^2) - (Ru)^2; sum((x - Cu2).^2) - (Ru2)^2], 'eq', []);

utol = 0.00;
divtol = 1e-4;

xu1 = ((Ru)^2 - sum((x - Cu).^2) + utol);
xu2 = ((Ru2)^2 - sum((x - Cu2).^2) + utol);
xu = xu1*xu2;

X1 = struct('ineq', [sum((x - Cu).^2) - (Ru)^2; sum((x - Cu2).^2) - (Ru2)^2], 'eq', []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data-driven settings 

eps = 2;    % offline noise
epsw = 2;   % online noise

eg = 2;     % system model... for data generation
Cons.df = 3;
Cons.dg = 0;
Cons.n = n;     % num of states
Cons.T = 80;    % num of samples

% collect data
generated_data = getDataRND(eg,eps,Cons.T);
[A,B,xi] = getCons(generated_data,Cons);

N = [A B; -A -B];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end

% find nontrivial constraints
e = [eps*ones(Cons.n*Cons.T,1)+xi; eps*ones(Cons.n*Cons.T,1)-xi];
[N, e] = nontrivial_constraints(N, e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%polynomials

f = [x(2); -x(1) + (1/3)*x(1)^3 - x(2)];
g = [0; 1];

Ulim = 5;
Rlim = 10;

if TWO_DISK
    [h, ch, mh] = polynomial(x, d-4);
    rho = h*xu1*xu2;
else
    [h, ch, mh] = polynomial(x, d-2);
    rho = h*xu1;
end


% [psi, cpsi, mpsi] = polynomial(x, d);
[psi, cpsi, mpsi] = polynomial(x, d, 1);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% when we have offline noise for robust data-driven

[~,~,vf] = polynomial(vars,Cons.df,1);
[~,~,vg] = polynomial(vars,Cons.dg,0);

dmax = max(d+Cons.dg-1, d+Cons.df-1);

Y = zeros(1, size(N,1), 'like', sdpvar);
for i = 1:size(N,1)
    Y(i) = polynomial(vars, dmax, 0);
end

D1 = [kron(jacobian(rho,vars),vf') kron(jacobian(psi,vars),vg')];
D2 = [rho*[jacobian(vf,x(1))' jacobian(vf,x(2))'] psi*[jacobian(vg,x(1))' jacobian(vg,x(2))']];

k = coefficients(-D1-D2-Y*N, vars);

% when we have online noise for robust control

YY = zeros(1,2*n,'like',sdpvar); % [Y, YY]
for j = 1:2*n
    YY(j) = polynomial(vars, dmax, 0);
end

kk = coefficients(-jacobian(rho,vars) -YY*[eye(n);-eye(n)], vars);
ee = [epsw*ones(n, 1); epsw*ones(n, 1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%psatz constraints

%radial limiter
limiter = (Ulim+(x'*x)/Rlim);
rhoscale = 0;
% limiter = Ulim;




[p_0, cons_0, coeff_0] = constraint_psatz(rho, X0, vars, d);
[p_u1, cons_u1, coeff_u1] = constraint_psatz(-rho, Xu1, vars, d);
[p_u2, cons_u2, coeff_u2] = constraint_psatz(-rho, Xu2, vars, d);
[p_top, cons_top, coeff_top] = constraint_psatz(limiter*rho - psi, X, vars, d);
[p_bot, cons_bot, coeff_bot] = constraint_psatz(limiter*rho + psi, X, vars, d);

coeff_Y = coefficients(Y,vars);
cons_Y = sos(Y);
coeff_YY = coefficients(YY,vars);
cons_YY = sos(YY);


if SEPARATE_YE % using Ye<0 and YYee<0:
    [p_div, cons_div, coeff_div] = constraint_psatz(-Y*e, X, vars, dmax);
    [p_w, cons_w, coeff_w] = constraint_psatz(-YY*ee, X, vars, dmax);

    coeff = [coeff_div; coeff_w; coeff_0; coeff_u1; coeff_u2; coeff_top; coeff_bot; ...
        ch; cpsi; coeff_Y; coeff_YY];
    cons = [cons_div;cons_w; cons_0; cons_u1; cons_u2; cons_top; cons_bot; ...
        k==0; kk==0; cons_Y; cons_YY];

else    % using Ye+YYee<0:
    [p_div, cons_div, coeff_div] = constraint_psatz(-Y*e-YY*ee, X, vars, dmax);

    coeff = [coeff_div; coeff_0; coeff_u1; coeff_u2; coeff_top; coeff_bot; ...
        ch; cpsi; coeff_Y; coeff_YY];
    cons = [cons_div; cons_0; cons_u1; cons_u2; cons_top; cons_bot; ...
        k==0; kk==0; cons_Y; cons_YY];
end

opts = sdpsettings('solver', 'mosek');

sol = solvesos(cons, 0, opts, coeff);


if sol.problem== 0
%% recovery
[crho, mrho] = coefficients(rho, vars);
rho_rec = value(crho)'*mrho;
psi_rec = value(cpsi)'*mpsi;
h_rec = value(ch)'*mh;

rhof = polyval_func(rho_rec, vars);
psif = polyval_func(psi_rec, vars);

boxlim = 1000;
ode_opts = odeset('RelTol', 1e-7, 'Events', @(t, x) box_event(t, x, boxlim));


uf = @(x) psif(x)./rhof(x);
Ff = polyval_func(f, vars);
Gf = polyval_func(g, vars);
hf = polyval_func(h_rec, vars);

Fclosed = @(x) Ff(x) + uf(x).*Gf(x);


%% sampling 

% Nsample = 10;
% sampler = @() ball_sample(1, 2)'*R0 + C0;
% traj = cell(Nsample, 1);
% Tmax = 10;
% for i = 1:Nsample
%     x0_curr = sampler();
%     traj{i} = ode23(@(t, x) Fclosed(x), [0, Tmax], x0_curr, ode_opts);
%     traj{i}.u = uf(traj{i}.y);
% end

%% plotting (nonrobust/ddc)

colors=linspecer(3);
blp = 3;

figure(12)
clf
theta = linspace(0, 2*pi, 150);
hold on
plot(C0(1) + R0*cos(theta), C0(2) + R0*sin(theta), 'k', 'LineWidth', 2)
plot(Cu(1) + Ru*cos(theta), Cu(2) + Ru*sin(theta), 'color', colors(2, :), 'LineWidth', 2)
plot(Cu2(1) + Ru2*cos(theta), Cu2(2) + Ru2*sin(theta), 'color', colors(2, :), 'LineWidth', 2)
% fimplicit(@(x1, x2) rhof([x1; x2]), [-1,1,-1,1]*blp, 'color', colors(1, :))
%fsurf(@(x1, x2) rhof([x1; x2]), [-1,1,-1,1]*blp)
xlabel('$x_1$', 'interpreter', 'latex')
xlabel('$x_2$', 'interpreter', 'latex')
title('Density Generated by Equality Constraint', 'fontsize', 14)

% for i = 1:Nsample
%     plot(traj{i}.y(1, :), traj{i}.y(2, :), 'color', colors(3, :))
% end


xlim(blp*[-1,1])
ylim(blp*[-1,1])
axis square

% syms z1 z2
% vec_rho = monomials([z1;z2], 0:d);
% RHO = value(crho)'*vec_rho;
f1 = fcontour(@(x, y) hf([x; y]),'c','LineWidth', 2);
Rutol = sqrt(Ru^2+utol);
Ru2tol = sqrt(Ru2^2+utol);
plot(Cu(1) + Rutol*cos(theta), Cu(2) + Rutol*sin(theta), 'color', 'b', 'LineWidth', 2)
plot(Cu2(1) + Ru2tol*cos(theta), Cu2(2) + Ru2tol*sin(theta), 'color', 'b', 'LineWidth', 2)

f1.LevelList = 0;





% figure(14)
% clf
% hold on
% for i = 1:Nsample
%     plot(traj{i}.x, traj{i}.u, 'color', colors(1, :))
% end
% xl = xlim;
% % plot(xl, Ulim*[1, 1], 'color', 'k', 'LineWidth', 2)
% % plot(xl, -Ulim*[1, 1], 'color', 'k', 'LineWidth', 2)
% plot(xl, 0*[1, 1], ':k', 'LineWidth', 2)
% 
% xlabel('$t$', 'interpreter', 'latex')
% ylabel('$u(t)$', 'interpreter', 'latex')
% title('Radially Bounded Input Actuation', 'fontsize', 14)


end