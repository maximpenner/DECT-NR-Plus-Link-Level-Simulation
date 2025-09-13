%
%   STF is mapped to transmit stream 0.
%
%   physical_resource_mapping_STF_cell is always a cell(1,3)
%
%       physical_resource_mapping_STF_cell(1,1) contains the subcarrier indices as described in 5.2.2
%
%       physical_resource_mapping_STF_cell(1,2) contains the OFDM symbol index where the STF is placed, always 0    
%
%       physical_resource_mapping_STF_cell(1,3) contains the values for each corresponding subcarriers
%
%   example size for Figure 4.5-2 a):
%
%        physical_resource_mapping_STF_cell =
%
%          1×3 cell array
%
%            {14×1 double}    {[0]}    {14×1 double}
%
%
function [physical_resource_mapping_STF_cell] = STF(numerology, k_b_OCC, N_eff_TX, b)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;
    
    % copy all relevant variables
    N_b_OCC = numerology.N_b_OCC;

    %% 5.2.2
    
    % first tupel
    i = 0 : 1 : N_b_OCC/8-1;
    k_i_0 = k_b_OCC(i*4 + MATLAB_INDEX_SHIFT);
    
    % second tupel
    i = N_b_OCC/8 : 1 : N_b_OCC/4-1;
    k_i_1 = k_b_OCC(N_b_OCC/2 + 3 + (i-N_b_OCC/8)*4 + MATLAB_INDEX_SHIFT);
    
    % concat and transpose
    k_i = [k_i_0, k_i_1];
    k_i = k_i';
    
    % we save data in cells
    physical_resource_mapping_STF_cell = cell(1,3);
    physical_resource_mapping_STF_cell(1,1) = {k_i};
    physical_resource_mapping_STF_cell(1,2) = {0};
    
    % base sequences
    % y_0_b1 = [-1i, -1i, -1, -1, 1i, -1i, 1i, ...
    %           -1i, -1i, 1i, 1i, -1, -1i, -1];
    % y_0_b2 = [1,-1,1,1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,1,-1,1,1,1,1,-1,1,-1];
    % y_0_b4 = [y_0_b2, y_0_b2];
    % y_0_b8 = [y_0_b4, y_0_b4];
    % y_0_b12 = [y_0_b4, y_0_b4, y_0_b4];
    % y_0_b16 = [y_0_b8, y_0_b8];

    % base sequences update from ETSI TS 103 636-3 V1.4.1 (2023-01)
    y_0_b1 = exp(1i*pi/4) * [1, -1,1,1, -1,1,1, -1,1,1,1, -1, -1, -1];
    y_0_b2 = exp(1i*pi/4) * [-1,1, -1,1,1, -1,1,1, -1,1,1,1, -1,1, -1, -1, -1,1, -1, -1, -1,1,1,1, -1, -1, -1, -1];
    y_0_b4 = exp(1i*pi/4) * [-1, -1, -1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1, -1, 1, -1, ...
        1, 1, 1, -1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1, -1, 1, -1];

    y_0_b4_r = (2*mod(1:numel(y_0_b4),2) - 1) .* fliplr(y_0_b4);
    y_0_b8 = [y_0_b4 y_0_b4_r];

    y_0_b8_r = (2*mod(1:numel(y_0_b8),2) - 1) .* fliplr(y_0_b8);
    y_0_b16 = [y_0_b8 y_0_b8_r];

    i = (0:1: (12*14-1));
    y_0_b12 = y_0_b16(i + 2*14 + MATLAB_INDEX_SHIFT);
    
    % synchronization training data symbols
    % NOTE: cyclic rotation depending on the number of effective antennas is already included
    i = 0 : 1 : N_b_OCC/4-1;
    idx = mod( i + 2*log2(N_eff_TX) , numerology.N_b_OCC/4 );
    switch b
        case 1
            y_STF_0_i = y_0_b1(idx + MATLAB_INDEX_SHIFT);
        case 2
            y_STF_0_i = y_0_b2(idx + MATLAB_INDEX_SHIFT);
        case 4
            y_STF_0_i = y_0_b4(idx + MATLAB_INDEX_SHIFT);
        case 8
            y_STF_0_i = y_0_b8(idx + MATLAB_INDEX_SHIFT);
        case 12
            y_STF_0_i = y_0_b12(idx + MATLAB_INDEX_SHIFT);
        case 16
            y_STF_0_i = y_0_b16(idx + MATLAB_INDEX_SHIFT);
    end
    
    %% create output variables (complex transpose .')
    physical_resource_mapping_STF_cell(1,3) = {y_STF_0_i.'};
end
