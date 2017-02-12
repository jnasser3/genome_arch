function local_hic_data = get_local_hic_data(varargin)
% Subsets a long hi-c profile to a local one

params = {'full_hic_data',...
    'position',...
    'resolution',...
    'width',...
    'kr_normalize',...
    'normalization_vector',...
    'peak_normalize'};
dflts = {[],...
    [],...
    [],...
    [],...
    true,...
    [],...
    false};
arg = parse_args(params,dflts,varargin{:});

position_rounded = round(arg.position/(arg.resolution*1000))*(arg.resolution*1000);
idx1 = arg.full_hic_data(:,1) == position_rounded;
idx2 = arg.full_hic_data(:,2) == position_rounded;

data1 = arg.full_hic_data(idx1,:);
data2 = arg.full_hic_data(idx2,:);

%% Normalize data
if arg.kr_normalize
    data1 = normalize_hic_data(data1,arg.normalization_vector,arg.resolution*1000);
    data2 = normalize_hic_data(data2,arg.normalization_vector,arg.resolution*1000);
end

temp1 = data1(:,[2 3]);
temp2 = data2(:,[1 3]);
centered_data = vertcat(temp1,temp2);

%% Prune data to within width of the gene start 
idx3 = (centered_data(:,1) > (arg.position - arg.width) & centered_data(:,1) < (arg.position + arg.width));
local_hic_data = sortrows(centered_data(idx3,:));

%% Normalize wrt peak height
%Divide hic contacts by signal at promoter
if arg.peak_normalize && size(local_hic_data,1) > 0
    diff_to_tss = abs(local_hic_data(:,1) - position_rounded);
    [~,idx] = min(diff_to_tss);
    peak_norm_factor = local_hic_data(idx(1),2);
    local_hic_data(:,2) = local_hic_data(:,2) / peak_norm_factor;
end

%% Hack to remove duplicate entry at 'position'
[~,idx] = duplicates(local_hic_data(:,1));
local_hic_data(idx(1),:) = [];

end

