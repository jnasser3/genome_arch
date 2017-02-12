%% Plot Hi-C (Un-normalized) contacts for a genomic region
resolution = 25; %in kb
chromosome = 18;
width = 2e6;
position = 65184000;%tss_table{353,3};%[38268656 75262618 128748315];
gene_name = {'DSEL'};%tss_table{353,1};
cell_lines = {'HUVEC','K562','IMR90','GM12878_primary'};

ds = temp_get_hic_data_multiple_lines('position',position,...
    'gene_name',gene_name,...
    'chromosome',chromosome,...
    'cell_lines',cell_lines,...
    'resolution',resolution,...
    'width',width);
temp_compare_hic_maps(ds)

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
