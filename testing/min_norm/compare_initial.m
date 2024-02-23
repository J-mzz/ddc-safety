mbc_nonauto_psatz
close all
%run mbc_nonauto.m first
%compare the rational trajectory with the closed loop trajectory
syms z1 z2
vars = [z1; z2];
% rho_rec = polyval_func(rho, vars);
%div(rho(f + g u))

rhof = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2));
rhog = jacobian(rho,vars)*g + rho*(jacobian(g(1),z1)+jacobian(g(2),z2));
psig = jacobian(psi,vars)*g + psi*(jacobian(g(1),z1)+jacobian(g(2),z2));

% DATA_rational = testSys(Region,f,g,u,1e-3);

th0 = matlabFunction(rhof);
th1 = matlabFunction(rhog);
th2 = matlabFunction(psig);

% v = monomials(vars,0:d_r/2);
% tol =  v'*1e-8*eye(length(v))*v;
% tol = matlabFunction(tol);
tol1 = 1e-5;


P = min_norm_mbc_optimizer(th0, th1, tol1);

u_out = P([0; -2.5]);
% u_min_norm = @(z1, z2) P([z1; z2]);
u_min_norm = @(z1, z2) u_min_norm_mbc(P, z1, z2, U);

XX = [];
U_min = [];
U_rat = [];

% for i = 1:length(DATA_rational.Traj)
for i = 10:10
    xx = DATA_rational.Traj{i}(:,1);
    % uu_min = U(xx(1), xx(2));
    % uu_rat = u_min_norm(xx(1), xx(2));

    uu_rat = U(xx(1), xx(2));
    uu_min = u_min_norm(xx(1), xx(2));

    XX = [XX,  xx];
    U_min = [U_min, uu_min];
    U_rat = [U_rat, uu_rat];
end

fprintf('umin: %0.4f, urat: %0.4f \n', uu_min, uu_rat);

%% diagnose issues

%the rational controller is satisfying the divergence constraint in terms
%of (rho, psi). But the divergence constraint is not satisfied when applied
%as rhof + u_rat*rhog

divx = DIV(xx(1), xx(2));
rhox = RHO(xx(1), xx(2));
psix = PSI(xx(1), xx(2));
ux = U(xx(1), xx(2));

rhofx = th0(xx(1), xx(2));
rhogx = th1(xx(1), xx(2));
divx_rat = rhofx + rhogx * uu_rat;
divx_min = rhofx + rhogx * uu_min;


% plot(1:30, U_min,'m' , 1:30, U_rat,'b')
% legend('min-norm','rational')

% warning off
% DATA_min_norm = testSys(Region,f,g,u_min_norm,0, 'm',1);
% warning on
% 
% xlim([-5,2])
% ylim([-5,2])
%{
%% compare the input actuations
figure(23)
% hold on
% for i = 1:length(DATA_min_norm.Input)
% % for i = 1:1
%     urat = DATA_rational.Input{i};
%     unorm = DATA_min_norm.Input{i};
%     subplot(2,1, 1)
%     plot(DATA_rational.Time{i}, urat, 'b');
% 
%     ylabel('u(t)')
%     title('Rational Controller', 'FontSize', 16)
%     subplot(2,1,2)
%     plot(DATA_min_norm.Time{i}, unorm, 'b');
%     ylabel('u(t)')
%     title('Min Norm Controller', 'fontsize', 16)
% xlabel('time')
% end

subplot(2,1,1)
for i = 1:length(DATA_rational.Input)
    urat = DATA_rational.Input{i};
    plot(DATA_rational.Time{i}, urat, 'b');
    hold on
end
ylabel('u(t)')
title('Rational Controller', 'FontSize', 16)

subplot(2,1,2) 
for i = 1:length(DATA_min_norm.Input)
    unorm = DATA_min_norm.Input{i};
    plot(DATA_min_norm.Time{i}, unorm, 'm');
    hold on
end

ylabel('u(t)')
title('Min Norm Controller', 'fontsize', 16)
xlabel('time')
%}

%% 
function u_out = u_min_norm_mbc(P, z1, z2, U)
    [u_out, errorcode]= P([z1; z2]);
    u_rational = U(z1, z2);
end