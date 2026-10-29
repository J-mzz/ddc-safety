function [A,B,xi] = getCons(X,Cons)
% GETCONS  Data matrices of the consistency set.
%   [A, B, xi] = getCons(X, Cons) returns A, B and xi such that every system
%   x_dot = f(x) + g(x)*u compatible with the data satisfies
%   |A*vec(f) + B*vec(g) - xi| <= eps, where f and g are expanded in
%   monomials of degree 1 to Cons.df and 0 to Cons.dg (f(0) = 0).
%
%   X     data matrix, one sample per row: [u, t, x1, x2, dx1, dx2]
%   Cons  problem dimensions (n, T, df, dg)
%
%   Note: the rows of X are shuffled with randperm (uses the global random
%   stream).

dmin = 1;   % lowest monomial degree in f (f(0) = 0)

df = Cons.df;
dg = Cons.dg;
n = Cons.n;
T = Cons.T;

lf = nchoosek(df+n, n) -dmin;         % length of monomial vector for f(x)
lg = nchoosek(dg+n, n);         % ... for g(x)

syms x1 x2
vars = [x1,x2];

p = monomials(vars, dmin:df);      % phi
g = monomials(vars, 0:dg);      % gamma

p_fun = matlabFunction(p);
g_fun = matlabFunction(g);

A = []; B = [];
p_data  = zeros(T, lf);
ug_data = zeros(T, lg);

ind = randperm(size(X,1),T);
X = X(ind,:);

for i = 1:T
    p_data(i,:) = feval(p_fun, X(i,3), X(i,4));
    
    if dg == 0
        ug_data(i,:) = X(i, 1);
    elseif dg > 0
        ug_data(i, :) = X(i, 1)*feval(g_fun, X(i, 3), X(i, 4));
    end
    
    A = [A; kron(eye(n),p_data(i,:))];
    B = [B; kron(eye(n),ug_data(i,:))];
end

dx = X(1:T,[5 6])';
xi = dx(:);

end