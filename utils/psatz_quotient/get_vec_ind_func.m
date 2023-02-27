function m = get_vec_ind_func(zz,z,index,mr)

% m = getVecInd(zz,z,index,mr)
% find the position of element in z corresponding the total dictionary zz
% and position index.
% mr denotes the number of 0-1 variables


if isempty(z)
    error('Element need to be located cannot be empty!');
else
    [~,nz] = size(zz);
    if size(z,2) ~= nz
        error('z is in wrong size');
    end
    n = size(z,1);
    m = zeros(1,n);
    z(:,1:mr) = z(:,1:mr)>=1;
    for i = 1:n
        equ = bsxfun(@eq,zz,z(i,:));
        ind = sum(equ,2)==nz;
        m(i) = index(ind);
    end
%     m = arrayfun(@(i) index(sum(bsxfun(@eq,zz,z(i,:)),2)==nz), 1:n);
end