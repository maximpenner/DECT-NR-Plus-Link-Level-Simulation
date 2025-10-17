function [samples_antenna] = ofdm_signal_generation_Cyclic_prefix_insertion(antenna_streams_mapped, ...
                                                                            k_b_OCC, ...
                                                                            N_PACKET_symb, ...
                                                                            N_TX, ...
                                                                            N_eff_TX, ...
                                                                            N_b_DFT, ...
                                                                            u, ...
                                                                            N_b_CP, ...
                                                                            oversampling)

    % when oversampling (os), our signal must remain in the center of the spectrum
    os_idx_start = (oversampling-1)*N_b_DFT/2 + 1;
    os_idx_end = os_idx_start + N_b_DFT - 1;
    
    % when oversampling, the variables N_b_DFT and N_b_CP have to be scaled accordingly
    N_b_DFT = oversampling * N_b_DFT;
    N_b_CP = oversampling * N_b_CP;

    % Matlab's 'comm.OFDMModulator' input must be 3d
    mat_antenna_streams_mapped_matlab_compatible = zeros(N_b_DFT, N_PACKET_symb, N_TX);
    for i=1:1:N_TX

        temp = cell2mat(antenna_streams_mapped(i));
        
        % TS in Figure 4.5-2 and 4.5-3 for an DFT size of 64:
        %
        %   The highest frequency (on the right in the spectrum) has an index 31 at the very top of the Figure.
        %   In the array 'temp' this corresponds to the index 1.
        %
        %   The dc in the OFDM spectrum has an index 0 in the center of the Figure.
        %   In the array 'temp' this corresponds to the index 32.
        %
        %   The lowest frequency (on the left in the spectrum) has an index -32 at the very bottom of the Figure.
        %   In the array 'temp' this corresponds to the index 64.
        %
        % Matlab's 'comm.OFDMModulator' assumes the following:
        %
        %   The highest frequency (on the right in the spectrum) must be at index 64 in 'temp'.
        %   The dc frequency must be at index 33 in 'temp'.
        %   The lowest frequency (on the left in the spectrum) must be at index 1 in 'temp'.
        %
        % Therefore we need to flip the spectrum upside-down.
        mat_antenna_streams_mapped_matlab_compatible(os_idx_start : os_idx_end,:,i) = flipud(temp);
    end    
    
    T_u_symb_in_samples = (N_b_DFT*9)/8;

    % modulate
    switch u
        case 1
            
            % cyclic prefix length from 6.3.5
            N_b_CP_vec = N_b_CP*ones(1,N_PACKET_symb);
            N_b_CP_vec(1) = (N_b_DFT*3)/4;

            STF_in_samples = (T_u_symb_in_samples*14)/9;

            % modulate into time domain
            hMod = comm.OFDMModulator('FFTLength',N_b_DFT, ...
                                        'NumGuardBandCarriers', [0; 0], ...
                                        'InsertDCNull', false, ...
                                        'PilotInputPort', false, ...
                                        'CyclicPrefixLength', N_b_CP_vec, ... 'Windowing', false,'WindowLength', 1, ...
                                        'NumSymbols', N_PACKET_symb, ...
                                        'NumTransmitAntennas', N_TX);
            samples_antenna = hMod(mat_antenna_streams_mapped_matlab_compatible);            
            
            % % compatible code, can be used to verify modulation
            % if 1==1
            %     samples_antenna_compare = zeros(size(samples_antenna));
            %     for i=1:1:N_TX
            % 
            %         % extract f-domain
            %         antenna_streams_mapped_i                                = zeros(N_b_DFT, N_PACKET_symb);
            %         antenna_streams_mapped_i(os_idx_start : os_idx_end, :)  = cell2mat(antenna_streams_mapped(i));
            %         antenna_streams_mapped_i_ud                             = flipud(antenna_streams_mapped_i);
            %         antenna_streams_mapped_i_s                              = fftshift(antenna_streams_mapped_i_ud, 1);
            %         antenna_streams_mapped_i_t                              = ifft(antenna_streams_mapped_i_s, N_b_DFT, 1);
            % 
            %         % create stf
            %         stf_t = antenna_streams_mapped_i_t(:,1);
            %         stf_t_cp = [stf_t(end - (N_b_DFT*3)/4 + 1 : end); stf_t];
            % 
            %         % create data
            %         data_t = antenna_streams_mapped_i_t(:,2:end);
            %         data_t_cp = [data_t(end- N_b_CP + 1: end,:); data_t];
            % 
            %         % concat
            %         samples_antenna_compare(:,i) = [stf_t_cp; reshape(data_t_cp, numel(data_t_cp), 1)];
            %     end
            % 
            %     % check if samples are the same
            %     diff = samples_antenna - samples_antenna_compare;
            %     if max(max(abs(diff))) > 1e-8
            %         error('Matlab and compatible code do not show the same result.');
            %     end
            % end
            
            % shorten GI for u=1, see Figure 5.1-1
            GI_length_in_samples_now = N_b_CP + N_b_DFT;
            GI_length_in_samples_short = (T_u_symb_in_samples*4)/9;
            GI_length_in_samples_remove = GI_length_in_samples_now - GI_length_in_samples_short;
            samples_antenna(end-GI_length_in_samples_remove+1:end,:) = [];
            
        case {2,4,8}
            
            % cyclic prefix length from 6.3.5.
            % Matlab can only assign cyclic prefixes up to a length of N_b_DFT, but we need 5/4*N_b_DFT.
            % We assign a prefix of 1/4*N_b_DFT to the first STF symbol and add a second STF symbol without a CP.
            N_b_CP_vec = N_b_CP*ones(1,N_PACKET_symb);
            N_b_CP_vec(1) = (N_b_DFT*1)/4;  % first symbol with 1/4*N_b_DFT CP
            N_b_CP_vec(2) = 0;              % additional STF symbol with no CP
            
            STF_in_samples = T_u_symb_in_samples*2;

            % remove one GI symbol at the end
            mat_antenna_streams_mapped_matlab_compatible(:,end,:) = [];
            
            % add one additional STF symbol
            mat_antenna_streams_mapped_matlab_compatible = [mat_antenna_streams_mapped_matlab_compatible(:,1,:), mat_antenna_streams_mapped_matlab_compatible];
            
            % modulate into time domain
            hMod = comm.OFDMModulator('FFTLength',N_b_DFT, ...
                                        'NumGuardBandCarriers', [0; 0], ...
                                        'InsertDCNull', false, ...
                                        'PilotInputPort', false, ...
                                        'CyclicPrefixLength', N_b_CP_vec, ... 'Windowing', false,'WindowLength', 1, ...
                                        'NumSymbols', N_PACKET_symb, ...
                                        'NumTransmitAntennas', N_TX);
            samples_antenna = hMod(mat_antenna_streams_mapped_matlab_compatible);

            % % compatible code, can be used to verify modulation
            % if 1==1
            %     samples_antenna_compare = zeros(size(samples_antenna));
            %     for i=1:1:N_TX
            % 
            %         % extract f-domain
            %         antenna_streams_mapped_i                                = zeros(N_b_DFT, N_PACKET_symb);
            %         antenna_streams_mapped_i(os_idx_start : os_idx_end, :)  = cell2mat(antenna_streams_mapped(i));
            %         antenna_streams_mapped_i_ud                             = flipud(antenna_streams_mapped_i);
            %         antenna_streams_mapped_i_s                              = fftshift(antenna_streams_mapped_i_ud, 1);
            %         antenna_streams_mapped_i_t                              = ifft(antenna_streams_mapped_i_s, N_b_DFT, 1);
            % 
            %         % create stf
            %         stf_t = antenna_streams_mapped_i_t(:,1);
            %         stf_t_cp = [stf_t(end - (N_b_DFT*1)/4 + 1 : end); stf_t; stf_t];
            % 
            %         % create data
            %         data_t = antenna_streams_mapped_i_t(:,2:end-1);
            %         data_t_cp = [data_t(end- N_b_CP + 1: end,:); data_t];
            % 
            %         % concat
            %         samples_antenna_compare(:,i) = [stf_t_cp; reshape(data_t_cp, numel(data_t_cp), 1)];
            %     end
            % 
            %     % check if samples are the same
            %     diff = samples_antenna - samples_antenna_compare;
            %     if max(max(abs(diff))) > 1e-8
            %         error('Matlab and compatible code do not show the same result.');
            %     end
            % end
    end
    
    % 6.3.5
    % power scaling is always a little troublesome:
    %
    %   Each sinuoid, z_k_l*exp(j2pi...), has a power of one, so the total power of the OFDM signal would be N_OCC_l.
    %   When using the ifft, however, matlab adds a factor 1/N_b_DFT to each sinuoid (https://de.mathworks.com/help/matlab/ref/ifft.html), so the total power is N_OCC_l/N_b_DFT^2.
    %   When adding a factor sqrt(N_b_DFT) to each sinuoid, the total power becomes N_OCC_l/N_b_DFT. So when all subcarriers were used, i.e. N_OCC_l=N_b_DFT, total power would be one.
    %
    %   To get an average power of one for the entire OFDM signal even when only a subset of subcarriers is used, we have to add a factor N_b_DFT/sqrt(N_OCC_l).
    %   The total power then becomes N_OCC_l/N_b_DFT^2 * (N_b_DFT/sqrt(N_OCC_l))^2 = 1.
    
    % scale all symbols but STF
    N_OCC_l = numel(k_b_OCC);
    samples_antenna(STF_in_samples+1:end,:) = N_b_DFT/sqrt(N_OCC_l)*samples_antenna(STF_in_samples+1:end,:);
    
    % scale STF
    N_OCC_0 = N_OCC_l/4;
    samples_antenna(1:STF_in_samples,:) = N_b_DFT/sqrt(N_OCC_0)*samples_antenna(1:STF_in_samples,:);
    
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % The STF in inserted only into the first transmit stream.
    % When using multiple antennas, we have to boost the STF so it has the same average power as the other OFDM symbols.
    %samples_antenna(1:STF_in_samples,:) = sqrt(N_eff_TX)*samples_antenna(1:STF_in_samples,:);
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!

    [n_samples, n_ant] = size(samples_antenna);

    assert(n_samples == N_PACKET_symb*(N_b_DFT*9)/8);
    assert(n_ant == N_TX);
end
