%run mbc_nonauto.m first
%compare the rational trajectory with the closed loop trajectory

vars = [z1; z2];
% rho_rec = polyval_func(rho, vars);
%div(rho(f + g u))

rhof = jacobian(rho,vars)*f + rho*(jacobian(f(1),z1)+jacobian(f(2),z2));
rhog = jacobian(rho,vars)*g + rho*(jacobian(g(1),z1)+jacobian(g(2),z2));

% DATA_rational = testSys(Region,f,g,u,1e-3);

th0 = matlabFunction(rhof);
th1 = matlabFunction(rhog);

u_min_norm = @(z1, z2) min_norm_scalar([z1, z2], th0, th1);

DATA_min_norm = testSys(Region,f,g,u_min_norm,1e-3, 'm');


%% compare the input actuations
figure(23)
hold on
for i = 1:length(DATA_min_norm.Input)
% for i = 1:1
    urat = DATA_rational.Input{i};
    unorm = DATA_min_norm.Input{i};
    subplot(2,1, 1)
    plot(DATA_rational.Time{i}, urat, 'b');

    ylabel('u(t)')
    title('Rational Controller', 'FontSize', 16)
    subplot(2,1,2)
    plot(DATA_min_norm.Time{i}, unorm, 'b');
    ylabel('u(t)')
    title('Min Norm Controller', 'fontsize', 16)
xlabel('time')
end