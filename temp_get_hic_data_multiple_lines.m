function full_data_struct = temp_get_hic_data_multiple_lines(position,chromosome,cell_lines,resolution,width)

full_data_struct = struct([]);
for ii = 1:length(cell_lines)
    %% Read data
    hic_data = importdata(sprintf('~/Documents/genome_arch/data/hic/%s/%dkb_resolution_intrachromosomal/chr%d/MAPQGE30/chr%d_%dkb.RAWobserved',cell_lines{ii},resolution,chromosome,chromosome,resolution),'\t');
    hic_norm_vector = importdata(sprintf('~/Documents/genome_arch/data/hic/%s/%dkb_resolution_intrachromosomal/chr%d/MAPQGE30/chr%d_%dkb.KRnorm',cell_lines{ii},resolution,chromosome,chromosome,resolution),'\t');

    cell_struct = struct([]);
    for jj = 1:length(position)
        local_hic_data = get_local_hic_data(hic_data,position(jj),resolution,width,true,hic_norm_vector);

        this_cell_struct = struct('cell_line',cell_lines(ii),...
            'gene_position',position(jj),...
            'genome_coordinates',local_hic_data(:,1),...
            'hic_contacts',local_hic_data(:,2));

        cell_struct = [cell_struct; this_cell_struct];
    end

data_struct = struct('cell_line',cell_lines(ii),'hic_data',cell_struct);
full_data_struct = [full_data_struct; data_struct];
end

%% Make into a 2-d struct array
full_data_struct = [full_data_struct.hic_data];
end

