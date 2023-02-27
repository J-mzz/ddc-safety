function [mom_out] = sym_moments(var_count, order)
%SYM_MOMENTS monomial powers and moments in the case of a symmetry
%Inputs:
%  var_count: an array of length 4:
%               [order-2 symmetric;
%                order-2 standard;
%                symmetric;
%                standard]
%
%           the symmetry is f(bs,  b, xs, x) = f(-bs, b, -xs, x) 
%           (binary symmetric,  binary, symmetric, standard)
%
% order:    degree of relaxation, monomial powers
%           moment matrix will have moments of up to degree d=2*order
%
%Outputs:   struct mom_out with fields:
%   monom_int:      monom_int: monomials that are on the interior of the
%                   moment matrices
%   index:          Don't know what this does
%   index_order2:   A logical array signaling which monomials will be
%                   reduced by order-2 reduction
%   M_even:         indices for the even moments (even total degree of
%                   symmetric components)
%   M_odd:          indices for the even moments (odd total degree of
%                   symmetric components)
%   N_int:          Number of moments on the interior of the moment matrix
%   N_even:         Number of even monomials generating M_even
%   N_odd:          Number of odd  monomials generating M_odd
%   order:          order of relaxation
%% preliminary counting and indexing

if length(var_count)==1;
    var_count = [0;0;0;var_count];
end
Nstd = var_count(3)+var_count(4);
Nbin = var_count(1)+var_count(2);

ind_sym = [(1:var_count(1)), var_count(1)+var_count(2) + (1:var_count(3))];

d = 2*order;

mom_out = struct;

%% Form monomial basis 

%find the monomial basis
[basis, ~] = momentPowers(Nbin, Nstd, order);

%split monomial basis into even and odd by symmetric variables
basis_sym = sum(basis(:, ind_sym), 2);
basis_mod = mod(basis_sym, 2);

basis_even = basis(basis_mod==0, :);
basis_odd  = basis(basis_mod==1, :);


%% Find monomial powers inside the matrix
m_even = size(basis_even, 1);
temp_even = bsxfun(@plus,kron(ones(m_even,1),basis_even),kron(basis_even,ones(m_even,1)));


m_odd= size(basis_odd, 1);
temp_odd = bsxfun(@plus,kron(ones(m_odd,1),basis_odd),kron(basis_odd,ones(m_odd,1)));

%monomials indexed by moment matrix
monom_int = unique([temp_even; temp_odd], 'rows');
index = 1:size(monom_int, 1);

index_order2 = any(monom_int(:, 1:sum(var_count(1:2))) == 2, 2);

M_even = reshape(get_vec_ind_func(monom_int,temp_even,index,0),m_even,[]);

if isempty(basis_odd)
    M_odd = [];
else
    M_odd  = reshape(get_vec_ind_func(monom_int,temp_odd, index,0),m_odd, []);
end


% ZZ_even = unique(temp_even, 'rows');
mom_out.monom_int = monom_int;
mom_out.index = index;
mom_out.index_order2 = index_order2; 
mom_out.M_even = M_even;
mom_out.M_odd = M_odd;

%number of monomials
mom_out.N_int = index(end);
mom_out.N_even = m_even;
mom_out.N_odd = m_odd;
mom_out.order = order;

end

