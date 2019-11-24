% [out_st,out_col] = CDF_linest(n1,n2)
% n1: number of styles in total
% n2: number of shape, maximum 12
% line_st = ['o','d','s','p','h','^','>','<','v','x','+','*'];

function [out_st,out_col] = CDF_linest(n1,n2)

    line_st = ['o','d','s','p','^','>','h','<','v','x','+','*'];
    line_st = line_st(1:n2);
    num_st  = numel(line_st);
    num_col = fix(n1/num_st)+1;
    col  = distinguishable_colors(num_col);
    col2 = lines(num_col);
    if num_col <= 7,
        col = col2(1:num_col,:);
    elseif num_col <= 14,
        col_temp = distinguishable_colors(20);
        col_temp = col_temp([5 7 10 11 12 17 20],:);
        col = [lines(7); col_temp];
    else
        col(1:6,:) = col2(1:6,:);
    end
    out_st = repmat(line_st,1,num_col);
    out_col = [];
    for i = 1:num_col
        out_col = [out_col; repmat(col(i,:),num_st,1)];
    end
    out_st = out_st(1:n1);
    out_col = out_col(1:n1,:);
end
    
% CDF_colorbar(col_temp,[1 2 3],[1 2 3],'',1,15)