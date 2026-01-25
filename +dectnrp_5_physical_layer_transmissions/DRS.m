%
%   DRS is mapped across all transmit streams, but the mapping is different for each stream.
%
%   physical_resource_mapping_DRS_cell is always a cell(N_TS,4)
%
%       N_TS is the number of transmit streams
%
%       physical_resource_mapping_DRS_cell(x,1) contains the subcarrier indices for each OFDM symbol as described in 5.2.3
%
%       physical_resource_mapping_DRS_cell(x,2) contains the OFDM symbol index where the DRS is placed  
%
%       physical_resource_mapping_DRS_cell(x,3) contains the values for each corresponding subcarriers
%                                               the values are the same for each OFDM symbol, but not for each transmit stream
%
%       physical_resource_mapping_DRS_cell(x,4) contains the linear indices of each subcarrier/OFDM-symbol pair
%
%
%   example size (N_TS = 8):
%
%        physical_resource_mapping_DRS_cell =
%
%          8×4 cell array
%
%            {14×10 double}    {1×10 double}    {14×1 double}    {14×10 double}
%            {14×10 double}    {1×10 double}    {14×1 double}    {14×10 double}
%            {14×10 double}    {1×10 double}    {14×1 double}    {14×10 double}
%            {14×10 double}    {1×10 double}    {14×1 double}    {14×10 double}
%            {14×10 double}    {1×10 double}    {14×1 double}    {14×10 double}
%            {14×10 double}    {1×10 double}    {14×1 double}    {14×10 double}
%            {14×10 double}    {1×10 double}    {14×1 double}    {14×10 double}
%            {14×10 double}    {1×10 double}    {14×1 double}    {14×10 double}
%
%
function [physical_resource_mapping_DRS_cell] = DRS(numerology, k_b_OCC, N_TS, N_eff_TX, N_PACKET_symb, b)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;
    
    % copy all relevant variables
    N_b_DFT = numerology.N_b_DFT;
    N_b_OCC = numerology.N_b_OCC;

    %% 5.2.3
    
    if N_eff_TX <= 2
        N_step = 5;
    elseif N_eff_TX >= 4
        N_step = 10;
    else
        error('N_eff_TX is %d, cannot determine N_step according to 5.2.3.', N_eff_TX);
    end
    
    i = 0 : 1 : N_b_OCC/4-1;

    nof_OFDM_symbols_carying_DRS = floor(N_PACKET_symb/N_step);

    % There's an error in the standard:
    % When using 4 TS and N_PACKET_symb is 15, 25, 35, 45, etc., the number of OFDM symbols carying DRS symbols is not floor(N_PACKET_symb/N_step).
    % It is one additional symbol, see figure 4.5-3 in part 3.
    if N_step == 10 && mod(N_PACKET_symb, 10) ~= 0
        
        assert(mod(N_PACKET_symb, 5) == 0);
        
        nof_OFDM_symbols_carying_DRS = nof_OFDM_symbols_carying_DRS + 1;
    end

    n = 0 : 1 : (nof_OFDM_symbols_carying_DRS-1);
    
    % we save data in cells
    physical_resource_mapping_DRS_cell = cell(N_TS,3);
    
    % iterate over each transmit stream
    for t = 0:1:N_TS-1
        i_column = i';
        k_i = k_b_OCC(i_column*4 + mod( t + mod(n,2)*2 ,4) + MATLAB_INDEX_SHIFT);
        l = 1 + floor(t/4) + n*N_step;
        
        % matlab related (if n ist a 1x1 vector, matlab create a row vector instead of a column vector)
        if isscalar(n)
            k_i = k_i';
        end

        physical_resource_mapping_DRS_cell(t + MATLAB_INDEX_SHIFT, 1) = {k_i};
        physical_resource_mapping_DRS_cell(t + MATLAB_INDEX_SHIFT, 2) = {l};

        %% linear indices
        rows = dectnrp_util.index_conversion_TS_matlab(N_b_DFT,k_i);
        cols = repmat(l + MATLAB_INDEX_SHIFT,size(rows,1),1);
        linear_indices_matlab = sub2ind([N_b_DFT (l(end) + MATLAB_INDEX_SHIFT)], rows, cols);
        physical_resource_mapping_DRS_cell(t + MATLAB_INDEX_SHIFT, 4) = {linear_indices_matlab};
    end
    
    % base sequences
    y_b1 = [1,1,1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,1, ...
            -1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1];
        
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % increase power when using more antennas
    %y_b1 = y_b1*sqrt(N_eff_TX);
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
        
    y_b2 = [y_b1,y_b1];
    y_b4 = [y_b2, y_b2];
    y_b8 = [y_b4, y_b4];
    y_b12 = [y_b4, y_b4, y_b4];
    y_b16 = [y_b8, y_b8];     
    
    % iterate over each transmit stream
    % signal transmitted in DRS subcarrier
    for t = 0:1:N_TS-1
        
        %i=0:1:N_b_DFT/4-1;         error in ETSI TS 103 636-3 V1.1.1.1 (2020-07)
        i=0:1:N_b_DFT/4 - 1 - 2*b;
        i_column = i';
        
        if t<4          % antenna indices 0,1,2,3, standard says t<=4, probably an error
            prefac = 1;
        elseif t>=4     % antenna indices 4,5,6,7, standard says t>4, probably an error
            prefac = -1;
        end
        
        switch b
            case 1
                y_DRS_i = y_b1(4*i_column + mod(t,4) + MATLAB_INDEX_SHIFT);
            case 2
                y_DRS_i = y_b2(4*i_column + mod(t,4) + MATLAB_INDEX_SHIFT);
            case 4
                y_DRS_i = y_b4(4*i_column + mod(t,4) + MATLAB_INDEX_SHIFT);
            case 8
                y_DRS_i = y_b8(4*i_column + mod(t,4) + MATLAB_INDEX_SHIFT);
            case 12
                y_DRS_i = y_b12(4*i_column + mod(t,4) + MATLAB_INDEX_SHIFT);
            case 16
                y_DRS_i = y_b16(4*i_column + mod(t,4) + MATLAB_INDEX_SHIFT);
        end
        
        y_DRS_i = prefac*y_DRS_i;
        
        y_DRS_i = y_DRS_i';
        
        physical_resource_mapping_DRS_cell(t + MATLAB_INDEX_SHIFT, 3) = {y_DRS_i};
    end
end
