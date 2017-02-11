function local_hic_data = get_local_hic_data(full_hic_data,position,resolution,width,kr_normalize,normalization_vector,peak_normalize)

    position_rounded = round(position/(resolution*1000))*(resolution*1000);
    idx1 = full_hic_data(:,1) == position_rounded;
    idx2 = full_hic_data(:,2) == position_rounded;

    %% Normalize data
    data1 = full_hic_data(idx1,:);
    data2 = full_hic_data(idx2,:);
    
    if kr_normalize
        data1 = normalize_hic_data(data1,normalization_vector,resolution*1000);
        data2 = normalize_hic_data(data2,normalization_vector,resolution*1000);
    end
    
    temp1 = data1(:,[2 3]);
    temp2 = data2(:,[1 3]);
    centered_data = vertcat(temp1,temp2);

    %% Prune data to within width of the gene start 
    idx3 = (centered_data(:,1) > (position - width) & centered_data(:,1) < (position + width));
    local_hic_data = sortrows(centered_data(idx3,:));
    
    %% Normalize wrt peak height
    %Divide hic contacts by signal at promoter
    if peak_normalize && size(local_hic_data,1) > 0
        diff_to_tss = abs(local_hic_data(:,1) - position_rounded);
        [~,idx] = min(diff_to_tss);
        %find(idx)
        %local_hic_data(idx,:)
        peak_norm_factor = local_hic_data(idx(1),2);
        local_hic_data(:,2) = local_hic_data(:,2) / peak_norm_factor;
    end


end

