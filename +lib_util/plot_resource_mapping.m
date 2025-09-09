function [] = plot_resource_mapping(numerology, N_PACKET_symb, N_TS, N_SS,...
                                    physical_resource_mapping_STF_cell,...
                                    physical_resource_mapping_DRS_cell,...
                                    physical_resource_mapping_PCC_cell,...
                                    physical_resource_mapping_PDC_cell)
                                
    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;
    
    % copy all relevant variables
    N_b_DFT = numerology.N_b_DFT;
    
    % create matrix
    [mat_STF_DRS_PCC_PDC, mat_STF_DRS_PCC_PDC_all_streams] = lib_util.matrix_STF_DRS_PCC_PDC(N_b_DFT, N_PACKET_symb, N_TS, N_SS,...
                                                                                                physical_resource_mapping_STF_cell,...
                                                                                                physical_resource_mapping_DRS_cell,...
                                                                                                physical_resource_mapping_PCC_cell,...
                                                                                                physical_resource_mapping_PDC_cell);

    %% STF
    figure(1)
    clf()
    imagesc(cell2mat(mat_STF_DRS_PCC_PDC(1)))
    title('STF in Transmit Stream 0');
    plot_setup(N_b_DFT);

    %% go over each stream for DRS
    for t=0:1:N_TS-1
        figure(2+t)
        clf()
        imagesc(cell2mat(mat_STF_DRS_PCC_PDC(2 + t)))
        title('DRS in some Transmit Stream');
        plot_setup(N_b_DFT);
    end

    %% PCC
    figure(1+N_TS+1)
    clf()
    imagesc(cell2mat(mat_STF_DRS_PCC_PDC(1+N_TS+1)));
    title('PCC in Spatial Stream 0');
    plot_setup(N_b_DFT);
    
    %% PDC
    for t=0:1:N_SS-1
        figure(1+N_TS+1+1+t)
        clf()
        imagesc(cell2mat(mat_STF_DRS_PCC_PDC(1+N_TS+1+1+t)))
        title('PDC in any Spatial Stream');
        plot_setup(N_b_DFT);
    end

    %% sanity check
    if 1+N_TS+1+N_SS ~= numel(mat_STF_DRS_PCC_PDC)
        error('Incorrect number of figures.');
    end

    %% all
    figure(1+N_TS+1+N_SS+1)
    clf()
    imagesc(mat_STF_DRS_PCC_PDC_all_streams);
    title('Resource Mapping');
    plot_setup(N_b_DFT);
end

function plot_setup(N_b_DFT)

    set(gca, 'YTick', 0:2:N_b_DFT, 'YTickLabel', N_b_DFT/2-(0:2:N_b_DFT))

    ylabel('Subcarrier Index');
    xlabel('OFDM symbol index');
    axis image
    axis ij
    colormap jet
    caxis([0 4])
end

