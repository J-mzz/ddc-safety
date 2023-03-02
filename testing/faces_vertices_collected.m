x = sdpvar(2, 1);
w = sdpvar(size(N, 2), 1);
C = randn(2, size(N, 2));
W = [N*w <= e; abs(x)<=2; uncertain(w); (C*w).*x <= 1];
cons = W;
sol = optimize(cons, sum(x));


