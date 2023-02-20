function DATA = testSys(Region,f,g,u,epsw)
%% plot closed loop phase protrait

r0 = Region.r0;
c0 = Region.c0;     % initial

n = 10; % # of samples
sample = (1:n)/n * 2*pi;
x = c0(1) + sqrt(r0) * sin(sample);
y = c0(2) + sqrt(r0) * cos(sample);


FF = matlabFunction(f+g*u);
F = matlabFunction(f);
U = matlabFunction(u);
% G = matlabFunction(g);
%-------------------------Trajectory: ode45

tspan = [0:.1:1.2]; % time span
curr_ode_options =  odeset('RelTol', 1e-7, ...
    'AbsTol', 1e-8, 'MaxStep', 1e-5);

Noise = @(t, x) epsw*(2*rand(2,1)-ones(2,1));
FFF = @(t,x) FF(x(1),x(2)) + Noise(x(1),x(2));

for i = 1:n
    x0 = [x(4);y(4)];% [x(9);y(9)];
    X = ode15s(FFF, tspan, x0, curr_ode_options);   % u \in [-1,1]
    
    DataX{i} = X.y(1,:);
    DataY{i} = X.y(2,:);
    
% 	plot(X.y(1,:), X.y(2,:),'c','LineWidth',2);
    plot(X.y(1,:), X.y(2,:),'c');
end
DATA.X = DataX;
DATA.Y = DataY;
%-------------------------Trajectory: streamline

% q = 5;     % bound
% k = 0.1;   % step size
% [XX,YY] = meshgrid(-q:k:q, -q:k:q);
% 
% for i = 1:size(XX,1)
%     for j = 1:size(XX,2)
% %         noise_value  = epsw*(2*rand(2,1)-ones(2,1));
% %         noise_weight = norm(epsw*(2*rand(2,1)-ones(2,1)),inf)/...
% %             U(XX(i,j),YY(i,j));
%         temp = FF(XX(i,j),YY(i,j)) + epsw*(2*rand(2,1)-ones(2,1));
%         UU(i,j) = temp(1);
%         VV(i,j) = temp(2);
%     end
% end
% % V = V + subs(u,{x1 x2}, {X Y});
% 
% % quiver(XX,YY,UU,VV, 'autoscalefactor', 10)
% h = streamline(XX,YY,UU,VV,x,y,k);
% 
% set(h,'Color','c')
% 
% axis(1*[-q q -q q])
% xlabel('x_1', 'fontsize', 13)
% ylabel('x_2', 'fontsize', 13)
% set(gca,'fontsize', 13)

%-------------------------

xlim([-5 5])
ylim([-5 5])
axis square




%% plot phase portrait

plotpp(@(t,x)  F(x(1),x(2)) + g * U(x(1),x(2))) 

end