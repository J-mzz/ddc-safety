function Cons = setCons(eg)
% SETCONS  Problem dimensions and polynomial degrees of Example 1.
%   Cons.n     number of states        Cons.T     number of samples
%   Cons.drho  degree of rho           Cons.dpsi  degree of psi
%   Cons.df    degree of f             Cons.dg    degree of g

switch eg
    case 2
        Cons.n 	  = 2;	% # of state
        Cons.T    = 80;	% # of sample
        Cons.drho = 4;	% degree of rho
        Cons.dpsi = 4;  % degree of psi
        Cons.df   = 3;  % degree of f
        Cons.dg   = 0;  % degree of g
end

end
