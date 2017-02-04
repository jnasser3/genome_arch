function temp_compare_hic_maps(ds)

%Plot

for ii = 1:size(ds,1)
    figure;
    for jj = 1:size(ds,2)
        plot(ds(ii,jj).genome_coordinates,ds(ii,jj).hic_contacts,...
            'DisplayName',ds(ii).cell_line)
        hold on
    end
end


end

