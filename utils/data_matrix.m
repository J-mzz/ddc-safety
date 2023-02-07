function [A,B,xi] = data_matrix(d_f,d_g,n,T)
% build according to eq.(14)
global x_glo
s_f = nchoosek(d_f+n, n)-1;       % size of monomial vector note that 1 doesn't occur in f(x) so -1
s_g = nchoosek(d_g+n, n);
x1 = sym('x1', 'real');           % X1 is used as data
x2 = sym('x2', 'real');
vars = [x1; x2];
phi = monomials(vars, 1:d_f);     % monomial vector of f
gamma = monomials(vars, 0:d_g);   % monomial vector of g
phi_fun = matlabFunction(phi);
gamma_fun = matlabFunction(gamma);

A = [];
B = [];
phi_data = zeros(T, s_f);
ugamma_data = zeros(T, s_g);
% 
sampling = linspace(0,5,T);
for i = 1:T
    [~,ind(i)] = min(abs(x_glo(:,2)-sampling(i)));
end
for i = 1:T
    phi_data(i, :) = feval(phi_fun, x_glo(ind(i), 3), x_glo(ind(i), 4));
    %     ugamma_data(i, :) = x_glo(i, 1)*feval(gamma_fun, x_glo(i, 3), x_glo(i, 4));    % if d_g > 0
    ugamma_data(i,:) = x_glo(ind(i), 1);    % if d_g = 0
    A = [A; kron(eye(n), phi_data(i,:))];
    B = [B; kron(eye(n), ugamma_data(i,:))];
end
x_dot = x_glo(ind,5:6)';
xi = x_dot(:);
% 
% for i = 1:T
%     phi_data(i, :) = feval(phi_fun, x_glo(i, 3), x_glo(i, 4));
%     %     ugamma_data(i, :) = x_glo(i, 1)*feval(gamma_fun, x_glo(i, 3), x_glo(i, 4));    % if d_g > 0
%     ugamma_data(i,:) = x_glo(i, 1);    % if d_g = 0
%     A = [A; kron(eye(n), phi_data(i,:))];
%     B = [B; kron(eye(n), ugamma_data(i,:))];
% end
% 
% x_dot = x_glo(1:T, 5:6)';
% xi = x_dot(:);                 % \xi
end