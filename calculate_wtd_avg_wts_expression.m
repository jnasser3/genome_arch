function wts = calculate_wtd_avg_wts_expression(expr_table,target_line,pred_line,gene_name)
%calculates the weights to be used in a weighted average computation given
%a target line and a set of predictor lines. Based on baseline expression.

gene_idx = strcmp(expr_table{:,'gene_name'},gene_name);
%target_expr = expr_table{gene_idx,lower(sprintf('%s_avg_FPKM',target_line))};
target_expr = expr_table{gene_idx,target_line};

%cellfun(@(s) sprintf('%s_avg_FPKM',s),pred_line,'Uni',false)

% predictor_idx = ismember(lower(expr_table.Properties.VariableNames),...
%     lower(cellfun(@(s) sprintf('%s_avg_FPKM',s),pred_line,'Uni',false)));

predictor_idx = ismember(expr_table.Properties.VariableNames,pred_line);

pred_expr = expr_table{gene_idx,predictor_idx};

wts_temp = 1 ./ abs(target_expr - pred_expr + eps);
wts_temp = wts_temp ./ sum(wts_temp);

wts = struct([]);

for ii = 1:length(pred_line)
    %this_line = sprintf('%s_avg_FPKM',pred_line{ii});
    this_line = pred_line{ii};
    idx = strcmpi(this_line,expr_table.Properties.VariableNames(predictor_idx));
    wts(ii).cell_line = pred_line{ii};
    wts(ii).wt = wts_temp(idx);
end


end

