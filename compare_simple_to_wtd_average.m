function [wts,simple_average,wtd_average] = compare_simple_to_wtd_average(ds,expr,target_line,pred_line,gene)
%Compare predicting target lines's HiC profile with a simple average of
%predictor lines to a wtd average based on expression of a gene

% rnaseq2hic_cell_map = containers.Map(...
%     {'Gm12878_avg_FPKM','Huvec_avg_FPKM','K562_avg_FPKM','Nhek_avg_FPKM','Nhlf_avg_FPKM'},...
%     {'GM12878_primary','HUVEC','K562','NHEK','IMR90'});

%% Compute averages
predictor_lines = setdiff({ds.cell_line},target_line);
target_idx = strcmp({ds.cell_line},target_line);
pred_idx = ismember({ds.cell_line},pred_line);
ds_pred = ds(pred_idx);
simple_average = calculate_average_hic_profile_1d(ds,...
    'cell_lines',predictor_lines);
wts = calculate_wtd_avg_wts_expression(expr,target_line,predictor_lines,gene);
wtd_average = calculate_average_hic_profile_1d(ds,...
    'wts',wts,...
    'cell_lines',predictor_lines);

%% Compute errors
simple_max_error = max(abs(ds(target_idx).hic_contacts - simple_average.hic_contacts));
simple_l1_error = sum(abs(ds(target_idx).hic_contacts - simple_average.hic_contacts));
wtd_max_error = max(abs(ds(target_idx).hic_contacts - wtd_average.hic_contacts));
wtd_l1_error = sum(abs(ds(target_idx).hic_contacts - wtd_average.hic_contacts));

%% Plot Results
figure;
plot(ds(target_idx).genome_coordinates,ds(target_idx).hic_contacts,'DisplayName',target_line{1})
hold on
plot(simple_average.genome_coordinates,simple_average.hic_contacts,'DisplayName','simple_average')
hold on
plot(wtd_average.genome_coordinates,wtd_average.hic_contacts,'DisplayName','wtd_average')

title_str = sprintf(['%s \n',...
    'simple max error = %.2f ',...
    'simple l1 error = %.2f \n',...
    'wtd max error = %.2f ',...
    'wtd l1 error = %.2f '],...
    gene,simple_max_error,simple_l1_error,wtd_max_error,wtd_l1_error);
title(title_str);
xlabel('Genomic Coordinates')
ylabel('HiC Contacts (peak normalized)')
l = legend;
set(l,'Interpreter','none')
legend show


end

