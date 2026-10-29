function [f,g] = getSystem(eg,vars)
% GETSYSTEM  True system x_dot = f(x) + g*u of Example 1 (flow system).
%   [f, g] = getSystem(eg, vars) evaluates f at vars (symbolic, sdpvar or
%   numeric). eg = 2 is the system of Example 1, from Rantzer's "On
%   analysis and synthesis of safe control laws".

x1 = vars(1);
x2 = vars(2);

switch eg
    case 2
        f = [x2; -x1 + 1/3*x1^3 - x2];
end

g = [0; 1];

end
