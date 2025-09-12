function [] = plot_resource_mapping_before_antenna(antenna_streams_mapped)
    
    N_b_DFT = size(cell2mat(antenna_streams_mapped), 1);

    total_sum = abs(cell2mat(antenna_streams_mapped(1)));
    total_sum_angle = abs(angle(cell2mat(antenna_streams_mapped(1))));
    for i=2:1:numel(antenna_streams_mapped)
        total_sum = total_sum + abs(cell2mat(antenna_streams_mapped(i)));
        total_sum_angle = total_sum_angle + abs(angle(cell2mat(antenna_streams_mapped(i))));
    end

    %% all
    figure()
    clf()
    subplot(2,1,1)
    imagesc(total_sum);
    title('Resource Mapping after antenna mapping, absolute');
    plot_setup(N_b_DFT);

    subplot(2,1,2)
    imagesc(total_sum_angle);
    title('Resource Mapping after antenna mapping, angle');
    plot_setup(N_b_DFT);
end

function plot_setup(N_b_DFT)

    set(gca, 'YTick', 0:2:N_b_DFT, 'YTickLabel', N_b_DFT/2-(0:2:N_b_DFT))

    ylabel('Subcarrier Index');
    xlabel('OFDM symbol index');
    axis image
    axis ij
    colormap jet
    %caxis([0 4])
end
