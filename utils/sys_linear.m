function dxdt = sys_linear(t,x,eps,l)
global x_glo
global index
if nargin == 3
    l = 1;
end
u = (2*rand(1,1) - 1)*l;
% x = x/5;
% u = - (5798836403786911*x(1))/3887665202774533 - (5595281382630794*x(2))/3887665202774533;
% u = (11.0*x(1) + 14.0*x(2))/(7.3e-8*x(1)^6 + 3.0e-7*x(1)^5*x(2) + 7.8e-6*x(1)^5 + 1.2e-6*x(1)^4*x(2)^2 + 4.7e-5*x(1)^4*x(2) + 4.9e-4*x(1)^4 + 1.8e-6*x(1)^3*x(2)^3 + 9.0e-5*x(1)^3*x(2)^2 + 1.1e-3*x(1)^3*x(2) - 6.9e-3*x(1)^3 + 1.5e-6*x(1)^2*x(2)^4 + 8.7e-5*x(1)^2*x(2)^3 + 1.2e-3*x(1)^2*x(2)^2 - 0.029*x(1)^2*x(2) - 0.15*x(1)^2 + 6.8e-7*x(1)*x(2)^5 + 4.3e-5*x(1)*x(2)^4 + 6.2e-4*x(1)*x(2)^3 - 0.027*x(1)*x(2)^2 - 0.1*x(1)*x(2) + 0.98*x(1) + 1.6e-7*x(2)^6 + 1.0e-5*x(2)^5 + 2.0e-4*x(2)^4 - 0.011*x(2)^3 - 0.07*x(2)^2 + 2.1*x(2) - 3.8);
x_glo(index, 1) = u(1);
x_glo(index, 2) = t;
x_glo(index, 3) = x(1);
x_glo(index, 4) = x(2);
noise = eps*(2*rand(2,1)-1)*l;
A = [0.5591    0.8962
    0.2408    0.5210]-0.5*eye(2);
% B = [0.6234    0.9294
%     0.0158    0.6909];
B = [0;1];
dxdt = A*x + B*u + noise;
dxdt = dxdt*l;
% dxdt = A*x + [0;1]*u + noise;
% dxdt = [x(2) + u + noise(1); -x(1)+1/3*x(1)^3-x(2) + noise(2)];
% dxdt = [2*x(1)^3 + x(1)^2*x(2) - 6*x(1)*x(2)^2 + 5*x(2)^3 + noise(1); u + noise(2)]; % example 2
% dxdt = [-6*x(1)*x(2)^2 - x(1)^2*x(2) + 2*x(2)^3 + noise(1); x(2)*u + noise(2)]; % example 3
% mu = 1;
% dxdt = [x(2) + noise(1) ; mu*(1-x(1)^2)*x(2) - x(1) + u + noise(2)];   % example 4
% dxdt = [1/(1-x(1)^2)*(u + x(1)) + noise(1); x(1) + noise(2)];  % example 5
% dxdt = [-x(1)+x(1)*x(2)+noise(1); -x(2)+noise(2)];  % example 6

x_glo(index, 5) = dxdt(1);
x_glo(index, 6) = dxdt(2);
% x_glo(index, 7) = u(2);
index = index + 1;
end