function attr_struct = gencode_attribute_to_struct(attr)
%Converts a string representing gencode attributes to a struct

temp = strsplit(attr,';');

attr_struct = struct;

for ii = 1:length(temp)
    this_temp = temp{ii};
    space1 = strfind(this_temp,' ');
    this_fname = this_temp(1:(space1-1));
    attr_struct.(this_fname) = strrep(this_temp((space1+1):end),'"','');
end

end