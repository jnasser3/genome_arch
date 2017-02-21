function gtf_out = annotate_gtf_with_attribute(gtf_in,attribute)
%Takes as input a gtf file
%Exatracts the fields in the cell array attribute to the main struct

% annotobj = GTFAnnotation(gtf_in);
% gtf = getData(annotobj);

gtf_out = rmfield(gtf_in,'Attributes');

for ii = 1:length(gtf_in)
    this_attr = gencode_attribute_to_struct(gtf_in(ii).Attributes);
    
    for jj = 1:length(attribute)
        gtf_out(ii).(attribute{jj}) = this_attr.(attribute{jj});
    end
    
end


end

