function [p_reduce] = yalmip_order2_reduce(p, var_reduce, g_substitute)
%YALMIP_ORDER2_REDUCE Reduction of a polynomial p with respect to the
%polynomial g = var_reduce^2 - g_substitute. The polynomial g is order-2
%with respect to var_reduce. This is for a single order-2 variable, and is
%very simple code.
%Examples: p mod x^2 - 1: yalmip_order2_reduce(p, x, 1)
%          p mod x^2 + y^2- 4: yalmip_order2_reduce(p, x, 4 - y^2)
%
% Inputs:
%       p:            polynomial to reduce
%       var_reduce:   which variable is order-2 and should be substituted
%       g_substitute: the remainder of var_reduce^2 after reduction

%break down the input polynomial with respect to var_reduce
if length(p) > 1
    p_reduce = p;
    for i = 1:length(p)
        p_reduce(i) = yalmip_order2_reduce(p(i), var_reduce, g_substitute);
    end
    return
end

%check to see if var_reduce is active
pow_basic = getexponentbase(p, var_reduce);
if ~any(pow_basic)
    %var_reduce not found in this expression
    p_reduce = p;
else
    q = 1;
    while any(q)
        [cc, mm] = coefficients(p, var_reduce);

        %perform the order-2 reduction
    %     mm_reduce = mm;
        pow_reduce =  getexponentbase(mm, var_reduce);

        %quotient and remainder
        r = mod(pow_reduce, 2);
        q = (pow_reduce - r)./2;
        mm_reduce = (g_substitute.^q) .* (var_reduce.^r);

        % for i = 1:length(mm)
        %     %quotient and remainder of each monomial of 
        %     r = mod(pow_reduce(i), 2);
        %     q = (pow_reduce(i)-r)/2;
        %     subs_mult = g_substitute^q;
        % 
        %     if r
        %         subs_mult = subs_mult * var_reduce;
        %     end
        %     mm_reduce(i) = subs_mult;
        % end

        %output reduced polynomial
        p_reduce = cc'*mm_reduce;
        p = p_reduce;
    end
end
end

