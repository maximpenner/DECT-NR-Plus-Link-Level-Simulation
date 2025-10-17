function [antenna_streams] = Beamforming(transmit_streams, N_TX, codebook_index)

    % remember N_eff_TX = N_TS
    N_TS = numel(transmit_streams);
    [N_b_DFT, N_PACKET_symb]= size(cell2mat(transmit_streams(1)));
    
    % create one big 3d matrix for input
    transmit_streams_mat = zeros(N_b_DFT, N_PACKET_symb, N_TS);
    for i=1:1:N_TS
        transmit_streams_mat(:,:,i) = cell2mat(transmit_streams(i));
    end
    
    % create 3d matrix for output
    antenna_streams_mat = zeros(N_b_DFT, N_PACKET_symb, N_TX);

    % extract correct precoding matrix
    W = dectnrp_6_generic_procedures.Beamforming_W(N_TS, N_TX, codebook_index);
    
    % this code is optimized for speed, see equivalent code below
    for i=1:1:N_TX
        %W_ = reshape(W(i,:),1,1,N_TX);  % ...oh Matlab
        W_ = reshape(W(i,:),1,1,[]);
        W_ = repmat(W_,N_b_DFT,N_PACKET_symb);
        antenna_streams_mat(:,:,i) = sum(W_ .* transmit_streams_mat,3);
    end
    
    % % compatible, but much slower code
    % if 1== 1
    %     antenna_streams_mat_check = zeros(N_b_DFT, N_PACKET_symb, N_TX);
    % 
    %     % go over each symbol and each subcarrier
    %     for i=1:1:N_PACKET_symb
    %         for j=1:1:N_b_DFT
    % 
    %             % 6.3.4
    %             y_N_TS = reshape(transmit_streams_mat(j,i,:), N_TS, 1);
    %             z_N_TX = W * y_N_TS;
    % 
    %             % write into output matrix
    %             antenna_streams_mat_check(j,i,:) = z_N_TX;
    %         end
    %     end
    % 
    %     diff = abs(antenna_streams_mat - antenna_streams_mat_check);
    %     if max(max(max(diff))) > 10e-6
    %         error('y_PcX not equal.');
    %     end
    % end
    
    % convert output to cells
    antenna_streams = cell(N_TX,1);
    for i=1:1:N_TX
        antenna_streams(i) = {antenna_streams_mat(:,:,i)};
    end
end
