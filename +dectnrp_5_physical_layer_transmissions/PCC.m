%
%   PCC is mapped to spatial stream 0.
%
%   physical_resource_mapping_PCC_cell is always a cell(1,x)
%
%       physical_resource_mapping_PCC_cell(1,1)
%       physical_resource_mapping_PCC_cell(2,1)
%       physical_resource_mapping_PCC_cell(3,1)
%       ...
%       physical_resource_mapping_PCC_cell(x-1,1) contain subcarrier indices for all OFDM symbols where the PCC is placed
%
%       physical_resource_mapping_PCC_cell(1,x) contains the corresponding OFDM symbol indices
%
%       physical_resource_mapping_PCC_cell(1,x) contains linear indices for fast demapping
%
%   example size for Figure 4.5-2 a):
%
%        physical_resource_mapping_PCC_cell =
%
%          1×3 cell array
%
%            {42×1 double}    {56×1 double}    {1×2 double}
%
%
function [physical_resource_mapping_PCC_cell] = PCC(numerology, ...
                                                    k_b_OCC, ...
                                                    N_TS, ...
                                                    N_PACKET_symb, ...
                                                    physical_resource_mapping_STF_cell, ...
                                                    physical_resource_mapping_DRS_cell)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;
    
    % copy all relevant variables
    N_b_DFT = numerology.N_b_DFT;    
    
    % to use the algorithm for PCC, we write all STF and DRS into one matrix of the size of the packet
    [~, mat_STF_DRS_all_streams] = dectnrp_5_physical_layer_transmissions.matrix_STF_DRS_PCC_PDC(N_b_DFT, N_PACKET_symb, N_TS, [], ...
                                                                                                 physical_resource_mapping_STF_cell, ...
                                                                                                 physical_resource_mapping_DRS_cell, ...
                                                                                                 [], ...
                                                                                                 []);


    %% 5.2.4
    
    l=1;
    N_PCC_subc = 98;
    N_unalloc_subc = N_PCC_subc;
    U = 0;
    
    % we save data in one cell
    %
    %   PCC in mapped into first spatial stream
    %
    %   number of rows: always 1
    %   first cell:     subcarrier indices for one OFDM symbol
    %   second cell:    subcarrier indices for one OFDM symbol
    %   third cell:     subcarrier indices for one OFDM symbol
    %   ...
    %   last cell:      indices of all aforementioned OFDM symbols (see OFDM_symbol_indices)
    %    
    physical_resource_mapping_PCC_cell = cell(0);
    OFDM_symbol_indices = [];
    
    % counter for found PDC subcarriers
    cnt_PCC_subc = 0;    
    
    % matlab has no do-while-loop, see break condition at the end
    while 1
        
        % find all indices of unused subcarriers in current OFDM-symbol
        [k_i_matlab, ~] = find(~mat_STF_DRS_all_streams(:,l + MATLAB_INDEX_SHIFT));
        
        % convert to notation of technical specification
        k_i = dectnrp_util.index_conversion_TS_matlab(N_b_DFT, k_i_matlab);
        
        % remove any subcarriers not within k_b_OCC, namely guards and DC
        k_i = intersect(k_i, k_b_OCC);
        
        U = numel(k_i);
        
        % no free subcarriers in this symbol, move to next OFDM symbol
        if U == 0
            l = l + 1;
            continue;
        end
        
        % sort from lowest subcarrier to highest
        k_i = sort(k_i);        
        
        % add to set of subcarriers
        k_PCC_l = k_i;
        
        % Step 5)
        R_PCC = 7;

        assert(mod(U,R_PCC)==0);

        % Step 6)
        C_PCC = U/R_PCC;

        % Step 7)
        R_PCC_x_C_PCC = reshape(k_PCC_l, C_PCC, R_PCC);
        R_PCC_x_C_PCC = R_PCC_x_C_PCC';

        % Step 8)
        k_PCC_l = reshape(R_PCC_x_C_PCC, numel(R_PCC_x_C_PCC), 1);
        
        % we either use as many subcarrieres as we can in this symbol (numel(k_PCC_l)) or as many as we still need (N_unalloc_subc)
        k_PCC_l = k_PCC_l(1 : min(numel(k_PCC_l), N_unalloc_subc));

        % sort from lowest subcarrier to highest
        k_PCC_l = sort(k_PCC_l);

        % append subcarriers for this OFDM symbol
        physical_resource_mapping_PCC_cell = [physical_resource_mapping_PCC_cell {k_PCC_l}];
        cnt_PCC_subc = cnt_PCC_subc + numel(k_PCC_l);
        
        % remember current OFDM symbol
        OFDM_symbol_indices = [OFDM_symbol_indices, l];
        
        % do while loop exit condition: we have found more subcarriers in this symbol (U) than we still need (N_unalloc_subc)
        if U >= N_unalloc_subc
            break;
        end
        
        % we have added some more subcarriers in this symbol, so we reduce the number of subcarriers required in the next OFDM symbol
        N_unalloc_subc = N_unalloc_subc - U;
        
        % move to next of OFDM symbol
        l = l + 1;
    end
    
    %% create a vector with linear indices
    linear_indices_matlab = zeros(cnt_PCC_subc, 1);
    idx = 1;
    for i=1:1:numel(OFDM_symbol_indices)
        
        k_i = cell2mat(physical_resource_mapping_PCC_cell(i));
        n_k_i = numel(k_i);
        
        % linear indices
        rows = dectnrp_util.index_conversion_TS_matlab(N_b_DFT,k_i);
        cols = repmat(OFDM_symbol_indices(i) + MATLAB_INDEX_SHIFT, numel(k_i), 1);
        li_matlab = sub2ind([N_b_DFT N_PACKET_symb], rows, cols);
        
        linear_indices_matlab(idx : idx + n_k_i - 1) = li_matlab;
        
        idx = idx + n_k_i;
    end  
    
    assert(idx-1 == cnt_PCC_subc && idx-1 == 98);
    
    % write result
    physical_resource_mapping_PCC_cell = [physical_resource_mapping_PCC_cell {OFDM_symbol_indices} {linear_indices_matlab}];    
end
