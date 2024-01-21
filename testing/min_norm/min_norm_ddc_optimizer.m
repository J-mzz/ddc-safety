function [P] = min_norm_ddc_optimizer(DR,D1,D2,N,e,epsw)
%get a yalmip optimizer for the data-driven (robust) min-norm controller 
u_min = sdpvar(1);
Y = sdpvar(1, size(N,1),'full');

n=2;
YY = sdpvar(1, 2*n,'full');

z = sdpvar(2, 1);

%parameters
% DR = sdpvar(size(DR0, 1), size(DR0, 2), 'full');
% D1 = sdpvar(size(D10, 1), size(D10, 2), 'full');
% D2 = sdpvar(size(D20, 1), size(D20, 2), 'full');


F = [
    [D1(z(1),z(2)), D2(z(1),z(2))*u_min] + Y*N == 0,...
    Y*e <= 0,...
    Y >= 0
    ];

if epsw ~= 0
    F = [F,...
        DR(z(1),z(2))+YY*[eye(n);-eye(n)] == 0,...
        YY*[epsw*ones(n, 1); epsw*ones(n, 1)] <= 0,...
        YY >= 0
        ];
end

options = sdpsettings('solver', 'mosek', 'verbose', 0);
% optimize(F, norm(u_min,2), options);


P = optimizer(F, norm(u_min, 2), options,{z}, u_min);

% u = value(u_min);

end
