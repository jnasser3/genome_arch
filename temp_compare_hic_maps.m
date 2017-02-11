function temp_compare_hic_maps(ds)

%Plot

for ii = 1:size(ds,1)
    f = figure;
    for jj = 1:size(ds,2)
        plot(ds(ii,jj).genome_coordinates,ds(ii,jj).hic_contacts,...
            'DisplayName',ds(ii,jj).cell_line)
        title(ds(ii,jj).gene_name)
        hold on
    end
    set(f,'Units','normalized')
    set(f,'OuterPosition',[0 0 1 1])
    legend show
    grid on
end


end

