function [ca,cc] = ddc_safety_b_red(A,B,xi,n,T,d,lambda,eps,K)
%data driven safety with a reduced number of faces

%% get size
d = num2cell(d);
[d_a,d_b,d_c,d_f,d_g] = deal(d{:});
%% Initilization
sdpvar x1 x2
vars = [x1 x2];
[a, ca, va] = polynomial(vars, d_a, 0);      % degree 0~d_a
% [~, ~, va] = polynomial(vars, d_a, 0);      % degree 0~d_a
[~, ~, vb] = polynomial(vars, d_b, 2);      % degree 2~d_b   use quadratic b(x)
[~, cc, vc] = polynomial(vars, d_c, 1);      % degree 1~d_c   omit constant control  change to 0~d_c for the case like xu
[~, ~, vf] = polynomial(vars, d_f, 1);      % degree 1~d_f   omit constant rate of change
[~, ~, vg] = polynomial(vars, d_g, 0);      % degree 0~d_g
[s1,c1] = polynomial(vars, d_a-2);
[s2,c2] = polynomial(vars, d_a-2);

% d_e = d_a;
% [e, ce, ve] = polynomial(vars, d_e, 0); 
cb = 1*[1;1;1];
% %% define safe/unsafe set
x0 = 0.16 - (x1-4)^2 - (x2-2)^2;
xu = 0.16 - (x1-2)^2 - (x2-2)^2;
%% Solve using YALMIPxs
dmax = max(d_b+d_c+d_g-1, d_b+d_a+d_f-1);
dmin = min(d_b+d_c+d_g-1, d_b+d_a+d_f-1);

[A_red, b_red] = nontrivial_constraints([A B; -A -B], [eps*ones(2*T, 1)+xi; eps*ones(2*T, 1)-xi]);
m_red = length(b_red);

Y = zeros(1, m_red, 'like', sdpvar);    % the polynomial matrix Y
for i = 1:m_red
%         [Y(i)] = polynomial(vars, dmax, dmin);   % low degree can be dropped
    [Y(i)] = polynomial(vars, dmax, 0);   % see (20) for size
end

% a = ca'*va;
b = cb'*vb;
af = ca'*(va*vf');       % ca^T m \phi^T
cg = cc'*(vc*vg');       % cc^T m \gamma^T

lieterm = af + cg;

% D1 = b* jacobian(lieterm, x

D1 = b*[ jacobian(af', x1)', jacobian(af', x2)', jacobian(cg', x1)', jacobian(cg', x2)'];
D2 = -lambda*[ jacobian(b, x1)'*af jacobian(b, x2)'*af, jacobian(b, x1)'*cg, jacobian(b, x2)'*cg];


D3 = K*a*b;
% D3 = K*e*b;
k = coefficients(-(D1+D2) - Y*A_red, vars);   % coefficients of polynomial vector kr-kl

% define constraints
F = [k == 0,...
    sos(-Y*b_red+D3),...
    sos(Y),...
    sos(a-0.01-s1*x0),...     % since rho = a/b where b > 0
    sos(-a-0.01-s2*xu),...
%     sos(e - a),...
%     sos(e + a),...
%     sum(ca) <= -0.1,...
    sos(s1),sos(s2)];

% sos solver
options = sdpsettings('solver', 'mosek', 'verbose', 0);
[sol,V,Q,res] = solvesos(F, [], options, [coefficients(Y,vars);c1;c2;ca;cc]);
sol
ca = value(ca);
cc = value(cc);
end