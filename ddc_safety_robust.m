function out = ddc_safety_robust(Cons,Region,eps,epsw,options)
% find the data-driven safety controller, given data, consistency set,
% (un)safe regions.

kkk = 1;

A   = Cons.A;
B   = Cons.B;
xi  = Cons.xi;
n   = Cons.n;
T   = Cons.T;
drho = Cons.drho;
dpsi = Cons.dpsi;
df  = Cons.df;
dg  = Cons.dg;
N = Cons.N;
e = Cons.e;

r0 = Region.r0;
ru = Region.ru;
c0 = Region.c0;
cu = Region.cu;

% initialization
sdpvar x1 x2
vars = [x1,x2];

[rho,cr] = polynomial(vars,drho,0);
[psi,cp] = polynomial(vars,dpsi,kkk);

[~,~,vf] = polynomial(vars,df,kkk);
[~,~,vg] = polynomial(vars,dg,0);

[s1,c1] = polynomial(vars,drho-2);
[s2,c2] = polynomial(vars,drho-2);

% slackness
[~,~,v] = polynomial(vars, drho/2, 0);
tol1 =  v'*1e-8*eye(length(v))*v;       % slackness, as mu in (20) 
tol2 =  v'*1e-8*eye(length(v))*v;       % slackness

% define initial and unsafe region
x0 = r0 - (x1-c0(1))^2 - (x2-c0(2))^2;	% rho > 0
xu = ru - (x1-cu(1))^2 - (x2-cu(2))^2;	% rho < 0

% solve with yalmip
dmax = max(dpsi+dg-kkk, drho+df-kkk);

Y = zeros(1, size(N,1), 'like', sdpvar);      % the polynomial matrix Y
for i = 1:size(N,1)
    Y(i) = polynomial(vars, dmax, 0);	% see (20) for size
end

YY = zeros(1,2*n,'like',sdpvar);        % extra part for noise
for j = 1:2*n
    YY(j) = polynomial(vars, dmax, 0);
end

% (15) in the note
D1 = [kron(jacobian(rho,vars),vf') kron(jacobian(psi,vars),vg')];
D2 = [rho*[jacobian(vf,x1)' jacobian(vf,x2)'] psi*[jacobian(vg,x1)' jacobian(vg,x2)']];

% k = coefficients(-D1-D2-Y*[A B; -A -B], vars);
k = coefficients(-D1-D2-Y*N, vars);

ubound = 1; % |psi/rho| <= ubound 

F = [k == 0,...
%     sos(-Y*[eps*ones(n*T, 1)+xi; eps*ones(n*T, 1)-xi] +tol1),...
    sos(-Y*e + 2*tol1),... %2e-6
%     sos(-Y*e - rho*xu +2e-6),... %2e-6; 
    (sos(Y)):'y',...
    ...
%     kk == 0,...
%     sos(-YY*[epsw*ones(n, 1); epsw*ones(n, 1)] +tol1),...
%     sos(YY),...
    ...
    sos(rho -s1*x0 -1*tol1),...	% since rho = a/b where b > 0
    sos(-rho-s2*xu -1*tol1),...
    ...
    sos(ubound*rho -psi),...    % bound input u
    sos(ubound*rho +psi),...
    ...
    sos(s1),sos(s2)];

if epsw ~= 0
    kk = coefficients(-jacobian(rho,vars) -YY*[eye(n);-eye(n)], vars);

    F = [F,...
        kk == 0,...
        sos(-YY*[epsw*ones(n, 1); epsw*ones(n, 1)]),... % +2*tol1),...
        sos(YY)];
end

sol = solvesos(F, [], options, [coefficients(Y,vars);c1;c2;cr;cp]);
out.cr = value(cr);
out.cp = value(cp);
out.sol = sol;
out.F = F;
out.Y = value(coefficients(Y,vars));
out.D1 = value(coefficients(D1,vars));
out.D2 = value(coefficients(D2,vars));

end
