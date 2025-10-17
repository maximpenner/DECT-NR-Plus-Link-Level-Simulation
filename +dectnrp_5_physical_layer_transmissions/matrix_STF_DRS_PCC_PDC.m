%
%   This function generates matrices, each of the size of a packet, i.e. N_b_DFT x N_PACKET_symb.
%
%   The matrices are filled with numbers, each number representing one channel:
%
%       0 for unoccupied
%       1 for STF
%       2 for DRS
%       3 for PCC
%       4*N_SS for PDC
%
%   Output:
%
%   mat_STF_DRS_PCC_PDC_all_streams contains all channels over all streams.
%
%   mat_STF_DRS_PCC_PDC is a cell(2+N_TS+N_SS,1):
%
%       mat_STF_DRS_PCC_PDC(1) STF
%
%       mat_STF_DRS_PCC_PDC(2) ... mat_STF_DRS_PCC_PDC(2 + N_TS - 1) DRS
%
%       mat_STF_DRS_PCC_PDC(2 + N_TS) PCC
%
%       mat_STF_DRS_PCC_PDC(2 + N_TS + 1) ... mat_STF_DRS_PCC_PDC(2 + N_TS + 1 + N_SS - 1) PDC
%
%   Arguments:
%
%       physical_resource_mapping_STF_cell required
%       physical_resource_mapping_DRS_cell required
%       physical_resource_mapping_PCC_cell can be left empty
%       physical_resource_mapping_PDC_cell can be left empty
%
function [mat_STF_DRS_PCC_PDC, mat_STF_DRS_PCC_PDC_all_streams] = matrix_STF_DRS_PCC_PDC(N_b_DFT, N_PACKET_symb, N_TS, N_SS, ...
                                                                                         physical_resource_mapping_STF_cell, ...
                                                                                         physical_resource_mapping_DRS_cell, ...
                                                                                         physical_resource_mapping_PCC_cell, ...
                                                                                         physical_resource_mapping_PDC_cell)
    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;
    
    assert(N_TS == size(physical_resource_mapping_DRS_cell,1));
   
    % types as integers
    type_STF = 1;
    type_DRS = 2;
    type_PCC = 3;
    type_PDC = 4;

    mat_STF_DRS_PCC_PDC = cell(0);

    %% STF

    mat_STF_in_transmit_stream_0 = zeros(N_b_DFT, N_PACKET_symb);

    % extract
    k_i = cell2mat(physical_resource_mapping_STF_cell(1,1));
    l = cell2mat(physical_resource_mapping_STF_cell(1,2));

    % convert
    k_i_matlab = dectnrp_util.index_conversion_TS_matlab(N_b_DFT, k_i);

    % write into matrix
    mat_STF_in_transmit_stream_0(k_i_matlab, l + MATLAB_INDEX_SHIFT) = type_STF;

    % append
    mat_STF_DRS_PCC_PDC = [mat_STF_DRS_PCC_PDC {mat_STF_in_transmit_stream_0}];

    %% DRS

    % go over each transmit stream
    for t=0:1:N_TS-1
        
        mat_DRS_single_transmit_stream = zeros(N_b_DFT, N_PACKET_symb); 
        
        % extract
        k_i_vec = cell2mat(physical_resource_mapping_DRS_cell(t + MATLAB_INDEX_SHIFT,1));
        l_vec = cell2mat(physical_resource_mapping_DRS_cell(t + MATLAB_INDEX_SHIFT,2));

        for q = 1:1:numel(l_vec)

            % extract single OFDM symbol
            l = l_vec(q);
            k_i = k_i_vec(:, q);
            
            % convert
            k_i_matlab = dectnrp_util.index_conversion_TS_matlab(N_b_DFT, k_i);
            
            % write into matrix
            mat_DRS_single_transmit_stream(k_i_matlab, l + MATLAB_INDEX_SHIFT) = type_DRS;
        end

        % append
        mat_STF_DRS_PCC_PDC = [mat_STF_DRS_PCC_PDC {mat_DRS_single_transmit_stream}];
    end

    %% PCC

    mat_PCC_in_spatial_stream_0 = zeros(N_b_DFT, N_PACKET_symb);

    % can be skipped
    if numel(physical_resource_mapping_PCC_cell) ~= 0

        % extract
        l_vec = cell2mat(physical_resource_mapping_PCC_cell(1,end-1));

        for q = 1:1:numel(l_vec)

            % extract single OFDM symbol
            l = l_vec(q);
            k_i = cell2mat(physical_resource_mapping_PCC_cell(1,q));

            % convert
            k_i_matlab = dectnrp_util.index_conversion_TS_matlab(N_b_DFT, k_i);

            % write into matrices
            mat_PCC_in_spatial_stream_0(k_i_matlab, l + MATLAB_INDEX_SHIFT) = type_PCC;
        end

        % append
        mat_STF_DRS_PCC_PDC = [mat_STF_DRS_PCC_PDC {mat_PCC_in_spatial_stream_0}];
    end

    %% PDC

    % can be skipped
    if numel(physical_resource_mapping_PDC_cell) ~= 0

        % go over each spatial stream
        for t=0:1:N_SS-1

            mat_PDC_single_transmit_stream = zeros(N_b_DFT, N_PACKET_symb);

            % extract
            l_vec = cell2mat(physical_resource_mapping_PDC_cell(1,end-1));

            for q = 1:1:numel(l_vec)

                % extract single OFDM symbol
                l = l_vec(q);
                k_i = cell2mat(physical_resource_mapping_PDC_cell(1,q));

                % convert
                k_i_matlab = dectnrp_util.index_conversion_TS_matlab(N_b_DFT, k_i);

                % write into matrices
                mat_PDC_single_transmit_stream(k_i_matlab, l + MATLAB_INDEX_SHIFT) = type_PDC;
            end

            % append
            mat_STF_DRS_PCC_PDC = [mat_STF_DRS_PCC_PDC {mat_PDC_single_transmit_stream}];
        end
    end

    %% entire matrix
    mat_STF_DRS_PCC_PDC_all_streams = zeros(N_b_DFT, N_PACKET_symb);
    for i=1:1:numel(mat_STF_DRS_PCC_PDC)
        mat_STF_DRS_PCC_PDC_all_streams = mat_STF_DRS_PCC_PDC_all_streams + cell2mat(mat_STF_DRS_PCC_PDC(i));
    end

    assert(isempty(setdiff(unique(mat_STF_DRS_PCC_PDC_all_streams),[0,1,2,3,N_SS*4])));
end
