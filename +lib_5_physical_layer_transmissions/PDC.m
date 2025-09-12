%
%   PDC is mapped across all spatial streams 0. The mapping is the same for each stream.
%
%   physical_resource_mapping_PDC_cell is always a cell(1,x)
%
%       physical_resource_mapping_PDC_cell(1,1)
%       physical_resource_mapping_PDC_cell(1,2)
%       physical_resource_mapping_PDC_cell(1,3)
%       ...
%       physical_resource_mapping_PDC_cell(1,x-2) contain subcarrier indices for all OFDM symbols where the PDC is placed
%
%       physical_resource_mapping_PDC_cell(1,x-1) contains the corresponding OFDM symbol indices
%
%       physical_resource_mapping_PDC_cell(1,x) contains linear indices for fast demapping
%
%   These indices are the same for each spatial stream.
%
%   example size fo Figure 4.5-2 b):
%
%        physical_resource_mapping_PDC_cell =
%
%          1×7 cell array
%
%            {56×1 double}    {56×1 double}    {56×1 double}    {42×1 double}    {56×1 double}    {56×1 double}    {1×6 double}    {322×1 double}
%
%
function [physical_resource_mapping_PDC_cell, N_PDC_re] = PDC(u, ...
                                                              numerology, ...
                                                              k_b_OCC, ...
                                                              N_PACKET_symb, ...
                                                              N_TS, ...
                                                              N_eff_TX, ...
                                                              physical_resource_mapping_STF_cell, ...
                                                              physical_resource_mapping_DRS_cell, ...
                                                              physical_resource_mapping_PCC_cell)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;
    
    % copy all relevant variables
    N_b_DFT = numerology.N_b_DFT;
    N_b_OCC = numerology.N_b_OCC;

    [~, mat_STF_DRS_PCC_all_streams] = lib_util.matrix_STF_DRS_PCC_PDC(N_b_DFT, ...
                                                                       N_PACKET_symb, ...
                                                                       N_TS, ...
                                                                       [], ...
                                                                       physical_resource_mapping_STF_cell, ...
                                                                       physical_resource_mapping_DRS_cell, ...
                                                                       physical_resource_mapping_PCC_cell, ...
                                                                       []);

    %% 5.2.5
    
    switch u
        case 1
            N_GI_plus_STF_symb = 2;
        case {2,4}
            N_GI_plus_STF_symb = 3;            
        case 8
            N_GI_plus_STF_symb = 4;            
    end
    
    N_DF_symb = N_PACKET_symb - N_GI_plus_STF_symb;
    
    if N_eff_TX <= 2
        N_step = 5;
    elseif N_eff_TX >= 4
        N_step = 10;
    else
        error('N_eff_TX is %d, cannot determine N_step according to 5.2.3.', N_eff_TX);
    end    

    nof_OFDM_symbols_carying_DRS = floor(N_PACKET_symb/N_step);

    % this is an error in the standard
    if N_step == 10 && mod(N_PACKET_symb, 10) ~= 0
        
        assert(mod(N_PACKET_symb, 5) == 0);
        
        nof_OFDM_symbols_carying_DRS = nof_OFDM_symbols_carying_DRS + 1;
    end

    N_DRS_re = N_eff_TX * N_b_OCC/4 * nof_OFDM_symbols_carying_DRS;
    
    % according to 5.2.4
    N_PCC_re = 98;
    
    N_PDC_re = N_DF_symb * N_b_OCC - N_DRS_re - N_PCC_re;

    %% extract subcarrier and symbol index for PDC
    
    % we save data in one cell
    %
    %   PDC in mapped into subcarriers which are not used in transmit streams by DRS or PCC in spatial stream 0
    %
    %   number of rows: always 1
    %   first cell:     subcarrier indices for one OFDM symbol
    %   second cell:    subcarrier indices for one OFDM symbol
    %   third cell:     subcarrier indices for one OFDM symbol
    %   ...
    %   last cell:      indices of all aforementioned OFDM symbols (see OFDM_symbol_indices)
    %    
    physical_resource_mapping_PDC_cell = cell(0);
    OFDM_symbol_indices = [];
    
    % first OFDM symbol at index 0 is STF, so we have to start at index 1
    l = 1;
    
    % counter for found PDC subcarriers
    N_PDC_re_cnt = 0;
    
    % loop over each OFDM sybol in Data Field (DF)
    for q = 1:1:N_DF_symb
        
        % find all indices of unused subcarriers in current OFDM-symbol
        [k_i_matlab, ~] = find(~mat_STF_DRS_PCC_all_streams(:,l + MATLAB_INDEX_SHIFT));
        
        % convert to notation of technical specification
        k_i = lib_util.index_conversion_TS_matlab(N_b_DFT, k_i_matlab);
        
        % remove any subcarriers not within k_b_OCC, namely guards and DC
        k_i = intersect(k_i, k_b_OCC);
        
        % no free subcarriers in this symbol, move to next OFDM symbol
        if numel(k_i)== 0
            l = l + 1;
            continue;
        end
        
        % sort from lowest subcarrier to highest
        k_i = sort(k_i);

        % append subcarriers for this OFDM symbol
        physical_resource_mapping_PDC_cell = [physical_resource_mapping_PDC_cell {k_i}];
        N_PDC_re_cnt = N_PDC_re_cnt + numel(k_i);
                
        % remember current OFDM symbol
        OFDM_symbol_indices = [OFDM_symbol_indices, l];

        % move to next of OFDM symbol
        l = l + 1;
    end   

    assert(N_PDC_re == N_PDC_re_cnt, "incorrect number of N_PDC_re");
    
    %% create a vector with linear indices
    linear_indices_matlab = zeros(N_PDC_re_cnt, 1);
    idx = 1;
    for i=1:1:numel(OFDM_symbol_indices)
        
        k_i = cell2mat(physical_resource_mapping_PDC_cell(i));
        n_k_i = numel(k_i);
        
        % linear indices
        rows = lib_util.index_conversion_TS_matlab(N_b_DFT,k_i);
        cols = repmat(OFDM_symbol_indices(i) + MATLAB_INDEX_SHIFT, numel(k_i), 1);
        li_matlab = sub2ind([N_b_DFT N_PACKET_symb], rows, cols);
        
        linear_indices_matlab(idx : idx + n_k_i - 1) = li_matlab;
        
        idx = idx + n_k_i;
    end

    assert(idx-1 == N_PDC_re_cnt, "incorrect number of linear indices");
    
    %% write result
    physical_resource_mapping_PDC_cell = [physical_resource_mapping_PDC_cell {OFDM_symbol_indices} {linear_indices_matlab}];
end
