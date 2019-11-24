% out_var = CDC_merge_group(in_var,NN)
% 
% Summarize the input according to catogories that is represented by NN elements
% Input should be a series of matrices in cell, with individual matrix
%   ending on the column of group indicators
% 
% NN is the number of columns used to indicate grouping

function out_var = CDC_merge_group(in_var,NN)

    kind_sum_temp = [];
    for i = 1:numel(in_var)
        if ~isempty(in_var{i}),
            temp = in_var{i}(:,end-NN+1:end);
            N(i) = size(in_var{i},2);
            kind_sum_temp = [kind_sum_temp; temp ones(size(temp,1),1)*i];
        else
            N(i) = NN + 1;
        end
    end
    N = N-NN;
    N_cum = cumsum([0 N]);

    [KIND_UNI,~,JJ] = unique(kind_sum_temp(:,1:NN),'rows');
    RECORD_OBS  = nan(size(KIND_UNI,1),nansum(N));
    for i = 1:numel(in_var)
        if ~isempty(in_var{i}),
            RECORD_OBS(JJ(kind_sum_temp(:,NN+1) == i),N_cum(i)+1:N_cum(i+1)) = in_var{i}(:,1:end-NN);
        end
    end

    out_var = [RECORD_OBS KIND_UNI];
end
