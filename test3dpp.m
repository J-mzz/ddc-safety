%% phase portrait

len = -2:0.8:2;

[X,Y,Z] = meshgrid(len,len,len);

% U = -Y.*Z + 1;
% V = X.*Z - Y;
% W = Z.^2 - Z.^3;

% U = -X + Y.^2;
% V = -Y + Z.^2;
% W = Z - X.^2;

U = -X;
V = -X -Y-Z-X.*Z;
W = X.*Y +Y;

f = quiver3(X,Y,Z,U,V,W,'Color',[0.5,0.5,0.5]);
hold on

%% Trajectories
% Ft_noise = @(t, y) [-y(1)+y(2)^2; -y(2)+y(3)^2; y(3)-y(1)^2];
Ft_noise = @(t, y) [-y(1); -y(1)-y(2)-y(3)-y(1)*y(3); y(1)*y(2)+y(2)];

y0 = [-1;-1;-1];
ode_options = odeset('RelTol', 1e-7,'AbsTol', 1e-8, 'MaxStep', 1e-6);
tspan = 0:10;

for i = 1:1
    [t,y] = ode15s(Ft_noise,tspan,y0,ode_options);
    plot3(y(:,1),y(:,2),y(:,3),'r')
end