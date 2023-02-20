function [f,g] = getSystem(eg,vars)

x1 = vars(1);
x2 = vars(2);

switch eg
    case 1
        A = [-0.6409    0.8962
            0.2408   -0.6790]; 
        f = A*vars; % random linear system
    case 2
        f = [x2; -x1 + 1/3*x1^3 - x2]; % ranzter's "On Analysis and Synthesis of Safe Control Laws"
    case 3
        A = [0 1; -1 -.2]; % Umich: Control Tutorials for MATLAB and Simulink
        f = A*vars;
    case 4
        f = [x1-x1*x2; -x2+x1*x2]; % predator prey model
    case 5
        f = [x2; -x1 - x1^3 - x2]; % nonlinear systems, 3rd edition, page 173
end

g = [0; 1];
% g = [0; x1];
 
end