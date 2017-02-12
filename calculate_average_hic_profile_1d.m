function avg_contact = calculate_average_hic_profile_1d(ds,varargin)
%Average together multiple hi-c profiles into an aggregate profile

params = {'cell_lines',...
    'wts'};
dflts = {'',...
    struct([])};
arg = parse_args(params,dflts,varargin{:});

%Subset ds to cell lines of interest
if ~isempty(arg.cell_lines)
    idx = ismember({ds.cell_lines},arg.cell_lines);
    ds = ds(idx);
end

%weights vector
if isempty(arg.wts)
    temp_wts = 1 / length(ds);
    wts = struct('cell_line',{ds.cell_line},...
        'weights',temp_wts);
end

%Get support of avg_profile as union of all provided supports
avg_support = [];
for ii = 1:length(ds)
    avg_support = union(avg_support,ds(ii).genome_coordinates);
end

%Map all provided profiles to bigger support
%Use linear interpolation
for ii = 1:length(ds)
    ds(ii).hic_contacts_interp = interp1(ds(ii).genome_coordinates,ds(ii).hic_contacts,avg_support);
end
    
%compute average
assert(isequal({wts.cell_line},{ds.cell_line}));
avg_profile = [wts.weights] * ([ds.hic_contacts_interp])';

avg_contact = struct('hic_contacts',avg_profile,...
    'genome_coordinates',avg_support');

end

