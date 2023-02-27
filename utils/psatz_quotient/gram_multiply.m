function [sigma_out] = gram_multiply(mom_in, monom,  Gram)
%GRAM_MULTIPLY perform sum(Gram .* monomials, 'all') = sigma(x)
%Input:
%   mom_in: moment matrix structure of indexing
%   monom:  monomials associated to the polynomial (possibly
%           quotient-reduced)
%   Gram:   PSD matrix variable from YALMIP

G_vec = accumarray(reshape(mom_in.M_even, [], 1), reshape(Gram, [], 1));

sigma_out = G_vec'*monom;

end

