function out_gtf = annotate_gtf_with_tss(gtf_in)
%Takes as input a gtf (Gencode) annotation struct.
%Outputs a gtf with relevant fields. 
%Calculates transcription start site based on strand

% annotobj = GTFAnnotation(gtf_in);
% gtf = getData(annotobj);

%tss_pos = cellfun( @get_tss, {gtf.Start}, {gtf.Stop}, {gtf.Strand});

out_gtf = gtf_in;
%[out_gtf(:).tss_pos] = tss_pos;
% out_gtf = arrayfun(@(t) setfield(t,'tss',tss_pos), gtf)
%out_gtf = arrayfun(@(t) setfield(t,'tss',get_tss(t.Start,t.Stop,t.Strand)), gtf);

for ii = 1:length(gtf_in)
    out_gtf(ii).tss = get_tss(gtf_in(ii).Start,gtf_in(ii).Stop,gtf_in(ii).Strand);
end

end


function tss = get_tss(start_pos,end_pos,strand)

if strcmp(strand,'+')
    tss = start_pos;
elseif strcmp(strand,'-')
    tss = end_pos;
else
    error('invalid strand')
end

end