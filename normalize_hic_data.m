function normalized_hic_data = normalize_hic_data(raw_hic_data,norm_vector,resolution)
%Normalizes Hi-C contact data via the norm_vector.
%
%As described in Rao et al 2014 (supplemental data README)

norm_idx = (raw_hic_data(:,[1 2]) / resolution + 1);
norm_factor1 = norm_vector(norm_idx(:,1));
norm_factor2 = norm_vector(norm_idx(:,2));

normalized_hic_data = raw_hic_data;
normalized_hic_data(:,3) = normalized_hic_data(:,3) ./ norm_factor1 ./ norm_factor2;
end

