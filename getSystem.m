function [f,g] = getSystem(eg,vars)

x1 = vars(1);
x2 = vars(2);

switch eg
    case 1
        A = [-0.6409    0.8962
            0.2408   -0.6790];
        f = A*vars;
    case 2
        f = [x2; -x1 + 1/3*x1^3 - x2];
    case 3
        A = [0 1; -5 -5];
        f = A*vars;
end

g = [0; 1];
 
end