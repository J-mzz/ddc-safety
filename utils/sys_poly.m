function dxdt = sys_poly(t,x,eps)
global x_glo
global index
u = 2*rand(1,1) - 1;
x_glo(index, 1) = u(1);
x_glo(index, 2) = t;
x_glo(index, 3) = x(1);
x_glo(index, 4) = x(2);
noise = eps*(2*rand(2,1)-1);
dxdt = [x(2) - x(1)^3 + x(1)^2 + noise(1); u + noise(2)];    % example1
% dxdt = [x(2)*q - x(1)^3*q^3 + x(1)^2*q^2 + noise(1); u + noise(2)];    % example1
% A = [0.5591    0.8962
%     0.2408    0.5210]-1*eye(2);
% B = [0.6234    0.9294
%     0.0158    0.6909];
% dxdt = A*x + B*u + noise;
% dxdt = A*x + [0 1]*u + noise;
% dxdt = [x(2) + u + noise(1); -x(1)+1/3*x(1)^3-x(2) + noise(2)];
% dxdt = [2*x(1)^3 + x(1)^2*x(2) - 6*x(1)*x(2)^2 + 5*x(2)^3 + noise(1); u + noise(2)]; % example 2
% dxdt = [-6*x(1)*x(2)^2 - x(1)^2*x(2) + 2*x(2)^3 + noise(1); x(2)*u + noise(2)]; % example 3
% mu = 1;
% dxdt = [x(2) + noise(1) ; mu*(1-x(1)^2)*x(2) - x(1) + u + noise(2)];   % example 4
% dxdt = [1/(1-x(1)^2)*(u + x(1)) + noise(1); x(1) + noise(2)];  % example 5
% dxdt = [-x(1)+x(1)*x(2)+noise(1); -x(2)+noise(2)];  % example 6

x_glo(index, 5) = dxdt(1);
x_glo(index, 6) = dxdt(2);
index = index + 1;
end