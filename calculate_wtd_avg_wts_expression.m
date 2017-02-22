function wts = calculate_wtd_avg_wts_expression(expr,target_line,pred_line,gene_name)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

gene_idx = strcmp(expr{:,'gene_name'},gene_name);
target_expr = expr{gene_idx,lower(sprintf('%s_avg_FPKM',target_line))};

cellfun(@(s) sprintf('%s_avg_FPKM',s),pred_line,'Uni',false)

predictor_idx = ismember(lower(expr.Properties.VariableNames),...
    lower(cellfun(@(s) sprintf('%s_avg_FPKM',s),pred_line,'Uni',false)));

pred_expr = expr{gene_idx,predictor_idx};

wts_temp = 1 ./ abs(target_expr - pred_expr);
wts_temp = wts_temp ./ sum(wts_temp);

wts = struct([]);

for ii = 1:length(pred_line)
    this_line = sprintf('%s_avg_FPKM',pred_line{ii});
    idx = strcmpi(this_line,expr.Properties.VariableNames(predictor_idx));
    wts(ii).cell_line = pred_line{ii};
    wts(ii).wt = wts_temp(idx);
    
end


end

