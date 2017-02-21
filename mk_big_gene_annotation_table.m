%% Inputs
rnaseq_lines = {'Gm12878','Huvec','K562','Nhek','Nhlf'};

%% Process Expression files
for ii = 1:length(rnaseq_lines)
    ii
    for jj = 1:2 %replicate number
        jj
        %~/Documents/genome_arch/data/rnaseq/caltech_ucsc_data
        %~/temp_small_caltech_rnaseq/
        filestr = fullfile('~/Documents/genome_arch/data/rnaseq/caltech_ucsc_data',...
            sprintf('wgEncodeCaltechRnaSeq%sR2x75Il200GeneGencV3cRep%dV3.gtf',rnaseq_lines{ii},jj));
        this_rnaseq_annot = GTFAnnotation(filestr);
        this_rnaseq_gtf = getData(this_rnaseq_annot);
        
        this_rnaseq_gtf = annotate_gtf_with_tss(this_rnaseq_gtf);
        this_rnaseq_gtf = annotate_gtf_with_attribute(this_rnaseq_gtf,...
                {'gene_name','FPKM','gene_type','gene_status'});
            
        % worst hack ever
        if jj == 1
            rep1 = this_rnaseq_gtf;
        elseif jj == 2
            rep2 = this_rnaseq_gtf;
        end

    end
    
    %Average FPKM over two replicates
    assert(isequal({rep1.Gene},{rep2.Gene}))
    avg = rep1;
    for kk = 1:length(avg)
        avg(kk).FPKM = (str2double(rep1(kk).FPKM) + str2double(rep2(kk).FPKM))/2;
    end
    
    %rename field
    fieldname = sprintf('%s_avg_FPKM',rnaseq_lines{ii});
    [avg.(fieldname)] = avg.FPKM;
    avg = rmfield(avg,'FPKM');
    
    if ii == 1
        T = struct2table(avg);
    else
        %make sure genes are in the same order
        assert(isequal(T{:,'Gene'},{avg.Gene}'),'Genes not in same order')
        
        temp = struct2table(keepfield(avg,fieldname));
        T = horzcat(T,temp);
    end
end