
function tss_table = temp_get_gene_and_tss
%% Read data
gene_pos = readtable('~/Documents/genome_arch/data/annot/gencode_genes_v19.txt','Delimiter','\t');
gene_annot = readtable('~/Documents/genome_arch/data/annot/gencode_attr_v19.txt','Delimiter','\t');

gene_join = join(gene_pos,gene_annot,...
    'RightKeys',{'transcriptId'},...
    'LeftKeys',{'name'});

%% Filter
chr_idx = strcmp(gene_join{:,'chrom'},'chr8');
protein_coding_idx = strcmp(gene_join{:,'transcriptType'},'protein_coding');
status_idx = strcmp(gene_join{:,'transcriptStatus'},'KNOWN');
all_idx = chr_idx & protein_coding_idx & status_idx;
this_table = gene_join(all_idx,:);

%% Make TSS Table
ugenes = unique(this_table{:,'geneName'});
tss_table1 = cell2table(cell(length(ugenes),2),...
    'VariableNames',{'gene_name','chr'});
tss_table2 = array2table(zeros(length(ugenes),1),...
    'VariableNames',{'approx_promoter_location'});
tss_table = horzcat(tss_table1,tss_table2);
% tss_table = struct2table(struct(length(ugenes),3),...
%     'VariableNames',{'gene_name','chr','approx_promoter_location'});


%all_idx = ones(1,length(ugenes));
for ii = 1:length(ugenes)
    this_idx = strcmp(this_table{:,'geneName'},ugenes{ii});
    
    if nnz(this_idx) > 1
        temp = find(this_idx);
        %all_idx(temp(2:end)) = 0;
        this_idx = temp(1);
    end
    
    if strcmp(this_table{this_idx,'strand'},'+')
        this_loc = this_table{this_idx,'txStart'};
    elseif strcmp(this_table{this_idx,'strand'},'-')
        this_loc = this_table{this_idx,'txEnd'};
    else
        error('strand must be either + or -')
    end
        
    tss_table{ii,'gene_name'} = this_table{this_idx,'geneName'};
    tss_table{ii,'chr'} = this_table{this_idx,'chrom'};
    tss_table{ii,'approx_promoter_location'} = this_loc;
end

end