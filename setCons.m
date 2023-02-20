function Cons = setCons(eg)

switch eg
    case {1,3,4}
        Cons.n 	  = 2;	% # of state
        Cons.T    = 40;	% # of sample
        Cons.drho = 4;	% degree of rho
        Cons.dpsi = 4;  % degree of psi
        Cons.df   = 1;  % degree of f
        Cons.dg   = 0;  % degree of g
    case {2,5}
        Cons.n 	  = 2;	% # of state
        Cons.T    = 80;	% # of sample
        Cons.drho = 4;	% degree of rho
        Cons.dpsi = 4;  % degree of psi
        Cons.df   = 3;  % degree of f
        Cons.dg   = 0;  % degree of g
end

end