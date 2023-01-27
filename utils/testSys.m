function testSys(Region,f,g,u)
% plot closed loop phase protrait


r0 = Region.r0;
c0 = Region.c0;     % initial

n = 20; % # of samples
sample = (1:n)/n * 2*pi;
x = c0(1) + sqrt(r0) * sin(sample);
y = c0(2) + sqrt(r0) * cos(sample);

tspan = [0,5]; % time span
eps = 1e-2; % noise level % 1e-2
curr_ode_options =  odeset('RelTol', 1e-7, ...
    'AbsTol', 1e-8, 'MaxStep', 0.01);

FF = matlabFunction(f+g*u+ eps*(2*rand(2,1)-1) );
FFF = @(t, x) FF(x(1), x(2));

for i = 1:n
    x0 = [x(i);y(i)];
    X = ode15s(FFF, tspan, x0, curr_ode_options);   % u \in [-1,1]

	plot(X.y(1,:), X.y(2,:),'c');
end

xlim([-5 5])
ylim([-5 5])
axis square


F = matlabFunction(f);
U = matlabFunction(u);

% plot phase portrait

% eps = 1e-2;
% RAN = @(t,x) eps*(2*rand(2,1)-1);

plotpp(@(t,x)  F(x(1),x(2)) + g * U(x(1),x(2))) 

end