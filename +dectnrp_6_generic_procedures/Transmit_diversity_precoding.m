function [y_PxC_ts] = Transmit_diversity_precoding(x_PxC_ss, N_TS)
    
    assert(ismember(N_TS, [2,4,8]));
    
    % always 98 elements in case of PCC
    x_PxC_ss = cell2mat(x_PxC_ss);
    n_x_PxC_ss = numel(x_PxC_ss);
    
    assert(mod(n_x_PxC_ss,2) == 0);
    
    % how many pairs of symbols do we have?, see 6.3.3.2
    i_max = n_x_PxC_ss/2;
    
    % modulo according to Table 6.3.3.2-1 to 6.3.3.2-3
    % switch N_TS
    %     case 2
    %         modulo = 1;
    %     case 4
    %         modulo = 6;
    %     case 8
    %         modulo = 12;
    % end
    
    % create vector of pairs (split real and imag) as shown in 6.3.3.2
    x_PxC_real = reshape(real(x_PxC_ss), 2, i_max);
    x_PxC_imag = reshape(imag(x_PxC_ss), 2, i_max);
    pair_of_symbols_vec = [x_PxC_real; x_PxC_imag];    
    
    % load precoding matrix in cells
    [Y, ~, ~] = dectnrp_6_generic_procedures.Transmit_diversity_precoding_Y(N_TS);

    [row,col,page] = size(Y);        

    % this block of code is optimized for fast calculation, for-loop below can be used instead
    pair_of_symbols_vec_serial = reshape(pair_of_symbols_vec,1,numel(pair_of_symbols_vec));
    pair_of_symbols_vec_serial = repmat(pair_of_symbols_vec_serial,row,1);
    Y_serial = reshape(Y, row, col*page);
    Y_serial = repmat(Y_serial, 1, ceil(numel(pair_of_symbols_vec)/(col*page)));
    Y_serial = Y_serial(:,1:size(pair_of_symbols_vec_serial,2));
    y_PxC_ts = pair_of_symbols_vec_serial.*Y_serial;
    y_PxC_ts = y_PxC_ts.';
    y_PxC_ts = reshape(y_PxC_ts,4,[]);
    y_PxC_ts = sum(y_PxC_ts,1);
    y_PxC_ts = reshape(y_PxC_ts,[],row);
    y_PxC_ts = y_PxC_ts.';

    % % compatible, but much slower code
    % if 1==1
    %     % preallocate output
    %     y_PxC_ts_check = zeros(2*N_TS,i_max);        
    % 
    %     % symbol by symbol
    %     for i=0:1:i_max-1
    % 
    %         % load correct Y matrix
    %         i_mod = mod(i,modulo);
    %         Y_mod = Y(:,:,i_mod+MATLAB_INDEX_SHIFT);
    % 
    %         % multiply
    %         i_matlab = i + MATLAB_INDEX_SHIFT;
    %         y_PxC_ts_check(:,i_matlab) = Y_mod*pair_of_symbols_vec(:,i_matlab);
    %     end
    % 
    %     diff = abs(y_PxC_ts - y_PxC_ts_check);
    %     if max(diff) > 10e-6
    %         error('y_PcX not equal.');
    %     end
    % end
    
    % rearrange output and convert to cell
    switch N_TS
        case 2
            y_PxC_TS_0 = y_PxC_ts([1,3],:);
            y_PxC_TS_0 = reshape(y_PxC_TS_0,numel(y_PxC_TS_0),1);
            
            y_PxC_TS_1 = y_PxC_ts([2,4],:);
            y_PxC_TS_1 = reshape(y_PxC_TS_1,numel(y_PxC_TS_1),1);
            
            y_PxC_ts = {y_PxC_TS_0, y_PxC_TS_1};
            
        case 4
            y_PxC_TS_0 = y_PxC_ts([1,5],:);
            y_PxC_TS_0 = reshape(y_PxC_TS_0,numel(y_PxC_TS_0),1);
            
            y_PxC_TS_1 = y_PxC_ts([2,6],:);
            y_PxC_TS_1 = reshape(y_PxC_TS_1,numel(y_PxC_TS_1),1);
            
            y_PxC_TS_2 = y_PxC_ts([3,7],:);
            y_PxC_TS_2 = reshape(y_PxC_TS_2,numel(y_PxC_TS_2),1);
            
            y_PxC_TS_3 = y_PxC_ts([4,8],:);
            y_PxC_TS_3 = reshape(y_PxC_TS_3,numel(y_PxC_TS_3),1);            
            
            y_PxC_ts = {y_PxC_TS_0; y_PxC_TS_1; y_PxC_TS_2; y_PxC_TS_3};
            
        case 8
            y_PxC_TS_0 = y_PxC_ts([1,9],:);
            y_PxC_TS_0 = reshape(y_PxC_TS_0,numel(y_PxC_TS_0),1);
            
            y_PxC_TS_1 = y_PxC_ts([2,10],:);
            y_PxC_TS_1 = reshape(y_PxC_TS_1,numel(y_PxC_TS_1),1);
            
            y_PxC_TS_2 = y_PxC_ts([3,11],:);
            y_PxC_TS_2 = reshape(y_PxC_TS_2,numel(y_PxC_TS_2),1);
            
            y_PxC_TS_3 = y_PxC_ts([4,12],:);
            y_PxC_TS_3 = reshape(y_PxC_TS_3,numel(y_PxC_TS_3),1);
            
            y_PxC_TS_4 = y_PxC_ts([5,13],:);
            y_PxC_TS_4 = reshape(y_PxC_TS_4,numel(y_PxC_TS_4),1);
            
            y_PxC_TS_5 = y_PxC_ts([6,14],:);
            y_PxC_TS_5 = reshape(y_PxC_TS_5,numel(y_PxC_TS_5),1);
            
            y_PxC_TS_6 = y_PxC_ts([7,15],:);
            y_PxC_TS_6 = reshape(y_PxC_TS_6,numel(y_PxC_TS_6),1);
            
            y_PxC_TS_7 = y_PxC_ts([8,16],:);
            y_PxC_TS_7 = reshape(y_PxC_TS_7,numel(y_PxC_TS_7),1);
            
            y_PxC_ts = {y_PxC_TS_0; y_PxC_TS_1; y_PxC_TS_2; y_PxC_TS_3; y_PxC_TS_4; y_PxC_TS_5; y_PxC_TS_6; y_PxC_TS_7};
    end 
end
