x = sdpvar(2, 1);
order = 2;
vars = x;
d = 2*order;
%start with the barrier

C0 = [1.5; 0];
R0 = 0.5;

Cu = [-1; -1];
Ru = 0.4;

X0 = struct('ineq', R0^2 - sum((x - C0).^2), 'eq', []);
Xu = struct('ineq', Ru^2 - sum((x - Cu).^2), 'eq', []);
X = struct('ineq', [], 'eq', []);
f = [x(2); -x(1) + (1/3)*x(1)^3 - x(2)];

[b, cb, mb] = polynomial(x, d);

Lb = jacobian(b, x)*f;

%psatz constraints

[p_lie, cons_lie, coeff_lie] = constraint_psatz(Lb, X, vars, d);
[p_0, cons_0, coeff_0] = constraint_psatz(b-1e-2, X0, vars, d);
[p_u, cons_u, coeff_u] = constraint_psatz(-b, Xu, vars, d);



coeff = [coeff_lie; coeff_0; coeff_u; cb];
cons = [cons_lie; cons_0; cons_u];
opts = sdpsettings('solver', 'mosek');
sol = solvesos(cons, 0, opts, coeff);

b_rec = value(cb)'*mb;
