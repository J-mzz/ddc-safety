function D = replacemoment(D, vars, Dict, m)
%% this code replace bilinear to linear, eg. m(2)*m(2) to m(4)
for k = 1:length(D)    % # of polynomials
    [c{k}, v{k}] = coefficients(D(k), vars);
    str = sdisplay(c{k});
    multi = regexp(str,'*');
    D_temp = 0;
    for i = 1:length(str)
        if isempty(multi{i})
            break;
        end
        for j = 1:length(multi{i})
            if multi{i}(j)-5 >= 1 & str{i}(multi{i}(j)-5) == 'm' & multi{i}(j)+5 <= length(str{i}) & str{i}(multi{i}(j)+5) == ')'
                ind1 = str2num(str{i}(multi{i}(j)+(-3:-2)));
                ind2 = str2num(str{i}(multi{i}(j)+(3:4)));
                IND = find(ismember(Dict, Dict(ind1,:)+Dict(ind2,:), 'rows'));
                old = str{i}(multi{i}(j)+(-5:5));
                new = char("m("+num2str(IND)+")");
            elseif multi{i}(j)-5 >= 1 & str{i}(multi{i}(j)-5) == 'm' & str{i}(multi{i}(j)+4) == ')'
                ind1 = str2num(str{i}(multi{i}(j)+(-3:-2)));
                ind2 = str2num(str{i}(multi{i}(j)+(3:3)));
                IND = find(ismember(Dict, Dict(ind1,:)+Dict(ind2,:), 'rows'));
                old = str{i}(multi{i}(j)+(-5:4));
                new = char("m("+num2str(IND)+")");
            elseif multi{i}(j)-4 >= 1 & str{i}(multi{i}(j)-4) == 'm' & multi{i}(j)+5 <= length(str{i}) & str{i}(multi{i}(j)+5) == ')'
                ind1 = str2num(str{i}(multi{i}(j)+(-2:-2)));
                ind2 = str2num(str{i}(multi{i}(j)+(3:4)));
                IND = find(ismember(Dict, Dict(ind1,:)+Dict(ind2,:), 'rows'));
                old = str{i}(multi{i}(j)+(-4:5));
                new = char("m("+num2str(IND)+")");
            elseif multi{i}(j)-4 >= 1 & str{i}(multi{i}(j)-4) == 'm' & str{i}(multi{i}(j)+4) == ')'
                ind1 = str2num(str{i}(multi{i}(j)+(-2:-2)));
                ind2 = str2num(str{i}(multi{i}(j)+(3:3)));
                IND = find(ismember(Dict, Dict(ind1,:)+Dict(ind2,:), 'rows'));
                old = str{i}(multi{i}(j)+(-4:4));
                new = char("m("+num2str(IND)+")");
            else
                continue;
            end
            OLD{j} = old;
            NEW{j} = new;
        end
        OLD = OLD(~cellfun('isempty',OLD));
        NEW = NEW(~cellfun('isempty',NEW));
        str{i} = replace(str{i},{OLD{:}},{NEW{:}});
        OLD = {};
        NEW = {};
        D_temp = D_temp + eval(str{i})*v{k}(i);
    end
    D(k) = D_temp;
end
end