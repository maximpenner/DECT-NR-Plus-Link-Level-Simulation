function [antenna_streams_mapped_rev, ...
          antenna_streams_mapped_rev_with_os_carriers] = ofdm_signal_generation_Cyclic_prefix_insertion_rev(samples_antenna_sto_cfo, ...
                                                                                                            k_b_OCC, ...
                                                                                                            N_PACKET_symb, ...
                                                                                                            N_RX, ...
                                                                                                            N_eff_TX, ...
                                                                                                            N_b_DFT, ...
                                                                                                            u, ...
                                                                                                            N_b_CP, ...
                                                                                                            oversampling)

    % remove the boost given to the packet either by beamforming or the scaling
    samples_antenna_sto_cfo = sqrt(N_eff_TX)*samples_antenna_sto_cfo;
                                                                                        
    % when oversampling, our signal must remain in the center of the spectrum
    os_idx_start = (oversampling-1)*N_b_DFT/2 + 1;
    os_idx_end = os_idx_start + N_b_DFT - 1;
    
    % when oversampling, the variables N_b_DFT and N_b_CP have to be scaled accordingly
    N_b_DFT = oversampling * N_b_DFT;
    N_b_CP = oversampling * N_b_CP;                                                                                        

    T_u_symb_in_samples = (N_b_DFT*9)/8;
    
    % switch to frequency domain
    switch u
        case 1
            
            % cyclic prefix length from 6.3.5
            N_b_CP_vec = N_b_CP*ones(1,N_PACKET_symb);
            N_b_CP_vec(1) = (N_b_DFT*3)/4;
            
            STF_in_samples = (T_u_symb_in_samples*14)/9;
            
            % 6.3.5
            % unscale all symbols but STF
            N_OCC_l = numel(k_b_OCC);
            samples_antenna_sto_cfo(STF_in_samples+1:end,:) = sqrt(N_OCC_l)/N_b_DFT*samples_antenna_sto_cfo(STF_in_samples+1:end,:);
            
            % unscale STF
            N_OCC_0 = N_OCC_l/4;
            samples_antenna_sto_cfo(1:STF_in_samples,:) = sqrt(N_OCC_0)/N_b_DFT*samples_antenna_sto_cfo(1:STF_in_samples,:);
            
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % we boosted the stf, remove this boost
            %samples_antenna_sto_cfo(1:STF_in_samples,:) = 1/sqrt(N_eff_TX)*samples_antenna_sto_cfo(1:STF_in_samples,:);
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            
            % we need to add some zeros to make it a multiple of N_b_DFT
            GI_length_in_samples = N_b_CP + N_b_DFT;
            GI_length_in_samples_short = (T_u_symb_in_samples*4)/9;
            GI_length_in_samples_add = GI_length_in_samples - GI_length_in_samples_short;
            samples_antenna_sto_cfo = [samples_antenna_sto_cfo; zeros(GI_length_in_samples_add, N_RX)];
            
            % OFDM demodulation, switch back to frequency domain
            demod = comm.OFDMDemodulator ('FFTLength',N_b_DFT, ...
                                          'NumGuardBandCarriers', [0; 0], ...
                                          'RemoveDCCarrier', false, ...
                                          'PilotOutputPort', false, ...
                                          'CyclicPrefixLength', N_b_CP_vec, ... 'Windowing', false,'WindowLength', 1, ...
                                          'NumSymbols', N_PACKET_symb, ...
                                          'NumReceiveAntennas', N_RX);
            antenna_streams_mapped_rev_matlab_compatible = step(demod,samples_antenna_sto_cfo);
            
        case {2,4,8}
            
            % cyclic prefix length from 6.3.5
            N_b_CP_vec = N_b_CP*ones(1,N_PACKET_symb);
            N_b_CP_vec(1) = (N_b_DFT*1)/4;
            N_b_CP_vec(2) = 0;
            
            STF_in_samples = T_u_symb_in_samples*2;
            
            % 6.3.5
            % unscale all symbols but STF
            N_OCC_l = numel(k_b_OCC);
            samples_antenna_sto_cfo(STF_in_samples+1:end,:) = sqrt(N_OCC_l)/N_b_DFT*samples_antenna_sto_cfo(STF_in_samples+1:end,:);
            
            % unscale STF
            N_OCC_0 = N_OCC_l/4;
            samples_antenna_sto_cfo(1:STF_in_samples,:) = sqrt(N_OCC_0)/N_b_DFT*samples_antenna_sto_cfo(1:STF_in_samples,:);
            
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % we boosted the stf, remove this boost
            %samples_antenna_sto_cfo(1:STF_in_samples,:) = 1/sqrt(N_eff_TX)*samples_antenna_sto_cfo(1:STF_in_samples,:);
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            % !!!!!!!!!!! NOT STANDARD COMPLIANT !!!!!!!!!!!
            
            % OFDM demodulation, switch back to frequency domain
            demod = comm.OFDMDemodulator ('FFTLength',N_b_DFT, ...
                                          'NumGuardBandCarriers', [0; 0], ...
                                          'RemoveDCCarrier', false, ...
                                          'PilotOutputPort', false, ...
                                          'CyclicPrefixLength', N_b_CP_vec, ... 'Windowing', false,'WindowLength', 1, ...
                                          'NumSymbols', N_PACKET_symb, ...
                                          'NumReceiveAntennas', N_RX);
            antenna_streams_mapped_rev_matlab_compatible = step(demod,samples_antenna_sto_cfo);
            
            % remove one preamble symbol which was added at the transmitter to fulfill the CP requirement
            antenna_streams_mapped_rev_matlab_compatible(:,1,:) = [];
            
            % add one additional zero GI symbol
            antenna_streams_mapped_rev_matlab_compatible = [antenna_streams_mapped_rev_matlab_compatible, zeros(N_b_DFT, 1, N_RX)];            
    end

    % convert to cell
    antenna_streams_mapped_rev_with_os_carriers = cell(N_RX,1);
    antenna_streams_mapped_rev                  = cell(N_RX,1);
    for i=1:1:N_RX
        antenna_streams_mapped_rev_with_os_carriers(i) = {flipud(antenna_streams_mapped_rev_matlab_compatible(1 : end,:,i))};
        antenna_streams_mapped_rev(i)                  = {flipud(antenna_streams_mapped_rev_matlab_compatible(os_idx_start : os_idx_end,:,i))};
    end
end
