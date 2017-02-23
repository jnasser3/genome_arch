function avg_contact = calculate_average_hic_profile_1d(ds,varargin)
%Average together multiple hi-c profiles into an aggregate profile

params = {'cell_lines',...
    'wts'};
dflts = {'',...
    struct([])};
arg = parse_args(params,dflts,varargin{:});

%Subset ds to cell lines of interest
if ~isempty(arg.cell_lines)
    idx = ismember({ds.cell_line},arg.cell_lines);
    ds_all = ds;
    ds = ds(idx);
end

%weights vector
if isempty(arg.wts)
    temp_wts = 1 / length(ds);
    wts = struct('cell_line',{ds.cell_line},...
        'wt',temp_wts);
else
    wts = arg.wts;
end

%Get support of avg_profile as union of all provided supports
%Need to use ds_all here. Need to include the target cell line in the union
%computation
avg_support = [];
for ii = 1:length(ds_all)
    avg_support = union(avg_support,ds_all(ii).genome_coordinates);
end

%Map all provided profiles to bigger support
%Use linear interpolation
for ii = 1:length(ds)
    ds(ii).hic_contacts_interp = interp1(ds(ii).genome_coordinates,...
        ds(ii).hic_contacts,...
        avg_support,...
        'linear',...
        'extrap');
end
    
%Put cell lines in same order
assert(isequal(sort({wts.cell_line}),sort({ds.cell_line})),'Cell line sets not the same!');
[~,loc] = ismember({ds.cell_line},{wts.cell_line});
wts = wts(loc);


assert(isequal({wts.cell_line},{ds.cell_line}),'Cell Lines not in same order!');
%compute average
assert(isequal({wts.cell_line},{ds.cell_line}),'Cell Lines not in same order!');
avg_profile = ([wts.wt] * ([ds.hic_contacts_interp])')';

avg_contact = struct('hic_contacts',avg_profile,...
    'genome_coordinates',avg_support',...
    'wts',wts);

end

