odefun = @(t, x) [x(2); -x(1) + 0.5*(1 - x(1)^2)*x(2)];
plotpp(odefun,'tspan', 30,...
                      'quivercolor', [0.6,0.6,0.6],...
                      'linecolor', [0.3,0.3,0.3])