%% Plot Hi-C (Un-normalized) contacts for a genomic region
resolution = 25; %in kb
width = 2e6;
gene_name = {'DSEL'};
cell_lines = {'HUVEC','K562','IMR90','GM12878_primary','NHEK'};

% Read rnaseq data and map cell line names
rnaseq2hic_cell_map = containers.Map(...
    {'Gm12878_avg_FPKM','Huvec_avg_FPKM','K562_avg_FPKM','Nhek_avg_FPKM','Nhlf_avg_FPKM'},...
    {'GM12878_primary','HUVEC','K562','NHEK','IMR90'});

expr = readtable('~/Documents/genome_arch/data/annotated_gene_expression.txt',...
    'Delimiter','\t');
for ii = 1:length(expr.Properties.VariableNames)
    if any(strcmp(expr.Properties.VariableNames(ii),rnaseq2hic_cell_map.keys))
        expr.Properties.VariableNames{ii} = rnaseq2hic_cell_map(expr.Properties.VariableNames{ii});
    end
end

idx = strcmp(gene_name,expr{:,'gene_name'});
position = expr{idx,'tss'};
chromosome = expr{idx,'Reference'}; chromosome = chromosome{1};

% position = 65184000;
% chromosome = 18;

ds = temp_get_hic_data_multiple_lines('position',position,... 
    'gene_name',gene_name,...
    'chromosome',chromosome,...
    'cell_lines',cell_lines,...
    'resolution',resolution,...
    'width',width);
%temp_compare_hic_maps(ds)

% Compare simple to wtd average
target_line = {'GM12878_primary'};
predictor_lines = setdiff(cell_lines,target_line);
[wts,simple_average,wtd_average] = compare_simple_to_wtd_average(ds,expr,target_line,predictor_lines,gene_name{1})

%% Compute averages
% 
% 
% target_idx = strcmp({ds.cell_line},target_line);
% ds_target = ds(target_idx);
% 
% ds_pred = ds(~target_idx);
% 
% %compute simple average profile
% simple_average = calculate_average_hic_profile_1d(ds_pred);
% 
% %compute wtd average based on expression
% %apply cell line map
% expr = readtable('~/Documents/genome_arch/data/annotated_gene_expression.txt',...
%     'Delimiter','\t');
% for ii = 1:length(expr.Properties.VariableNames)
%     if any(strcmp(expr.Properties.VariableNames(ii),rnaseq2hic_cell_map.keys))
%         expr.Properties.VariableNames{ii} = rnaseq2hic_cell_map(expr.Properties.VariableNames{ii});
%     end
% end
% %expr = table2struct(expr);
% 
% %'GM12878'
% %{'HUVEC','NHLF','K562','NHEK'}
% wts = calculate_wtd_avg_wts_expression(expr,target_line,predictor_lines,'DSEL');
% wtd_average = calculate_average_hic_profile_1d(ds_pred,'wts',wts);
% 
% temp_compare_hic_maps(ds)
% hold on
% plot(simple_average.genome_coordinates,simple_average.hic_contacts,'DisplayName','simple_average')
% hold on
% plot(wtd_average.genome_coordinates,wtd_average.hic_contacts,'DisplayName','wtd_average')
%% Read GTF file
% annotobj = GTFAnnotation('~/Documents/genome_arch/data/rnaseq/caltech_ucsc_data/wgEncodeCaltechRnaSeqGm12878R2x75Il200GeneGencV3cRep1V3.gtf');
% gtf = getData(annotobj);

%% Read in the file
% hic_data = importdata(sprintf('~/Documents/genome_arch/data/hic/%s/%dkb_resolution_intrachromosomal/chr%d/MAPQGE30/chr%d_%dkb.RAWobserved',cell_line,resolution,chromosome,chromosome,resolution),'\t');
% hic_norm_vector = importdata(sprintf('~/Documents/genome_arch/data/hic/%s/%dkb_resolution_intrachromosomal/chr%d/MAPQGE30/chr%d_%dkb.KRnorm',cell_line,resolution,chromosome,chromosome,resolution),'\t');
% 
% %% find coordinates
% position_rounded = round(position/(resolution*1000))*(resolution*1000);
% idx1 = hic_data(:,1) == position_rounded;
% idx2 = hic_data(:,2) == position_rounded;
% 
% %% Normalize data
% raw_data1 = hic_data(idx1,:);
% raw_data2 = hic_data(idx2,:);
% normalized_data1 = normalize_hic_data(raw_data1,hic_norm_vector,resolution*1000);
% normalized_data2 = normalize_hic_data(raw_data2,hic_norm_vector,resolution*1000);
% 
% temp1 = normalized_data1(:,[2 3]);
% temp2 = normalized_data2(:,[1 3]);
% centered_data = vertcat(temp1,temp2);
% 
% %% Prune data to within width of the gene start 
% idx3 = (centered_data(:,1) > (position - width) & centered_data(:,1) < (position + width));
% centered_data_pruned = sortrows(centered_data(idx3,:));
% %centered_data_pruned = centered_data(idx3,:);
% 
% %% Expand data to include regions with zero counts
% %full_genomic_loci = (position_rounded-width):resolution:(position_rounded+width);
% %missing_loci = setdiff(full_genomic_loci,centered_data_pruned(:,1))
% 
% %% Plot
% figure
% plot(centered_data_pruned(:,1),centered_data_pruned(:,2))
