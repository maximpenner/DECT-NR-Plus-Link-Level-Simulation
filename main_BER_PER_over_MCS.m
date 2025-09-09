clear all;
close all;

% This script calculates BERs and PERs for DECT-2020 New Radio packets over various wireless channels (AWGN, Rayleigh, Rician) for different MCSs.
% 
% For each MCS and each SNR, the same number of packets is calculated.
% Results are saved in the folder results/.
% When this script is finished, the scripts main_BER_PER_over_MCS_plot_PCC.m and main_BER_PER_over_MCS_plot_PDC.m can be used to plot the results.
%
% Executing this script as is should take only a few seconds to minutes, depending of the system and multi-core capabilities (parfor).

rng('shuffle');
%rng(1140598280);

warning('off');

if exist('results', 'dir')
    lib_util.clear_directory('results');
else
    mkdir('results');
end

fprintf('Starting at %s\n', datestr(now,'HH:MM:SS'));

% choose mcs to simulate and maximum number of harq retransmissions
mcs_index_vec = [1,2,3,4];      % part 3, Table A-1
max_harq_retransmissions = 0;

% simulation range for link level simulation
snr_db_vec_global = -10 : 1.0 : 30;
snr_db_vec_global = repmat(snr_db_vec_global,numel(mcs_index_vec),1);

% Packets per mcs and snr. Increase this number to get smoother curves.
n_packets_per_snr = 0.5e3;

% result container for PCC
n_bits_PCC_sent_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));        % BER uncoded
n_bits_PCC_error_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));       % BER uncoded
n_packets_PCC_sent_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));     % PER
n_packets_PCC_error_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));    % PER

% result container for PDC
n_bits_PDC_sent_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));        % BER uncoded
n_bits_PDC_error_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));       % BER uncoded
n_packets_PDC_sent_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));     % PER
n_packets_PDC_error_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));    % PER

bps_global = zeros(numel(mcs_index_vec), 1);        % bits per symbol, PDC only
tbs_global = zeros(numel(mcs_index_vec), 1);        % transport block size, PDC only

cnt = 1;
for mcs_index = mcs_index_vec

    % these variables need to be set before creating tx and rx
    mac_meta_tx.u = 1;                                  % mu = 1, 2, 4 or 8
    mac_meta_tx.b = 1;                                  % beta = 1, 2, 4, 8, 12 or 16
    mac_meta_tx.PacketLengthType = 0;                   % 0 for subslots, 1 for slots
    mac_meta_tx.PacketLength = 2;                       % min is 1, max is 16 according to Table 6.2.1-2a in part 4
    mac_meta_tx.tm_mode_0_to_11 = 0;                    % Table 7.2-1, mode determines wether transmission is closed loop or not, values range from 0 to 11
    mac_meta_tx.mcs_index = mcs_index;                  % Table A-1 in part 3, values range from 0 to 11
    mac_meta_tx.Z = 6144;                               % 5.3 -> so far only Z=6144 fully supported, 2048 only at TX, RX missing (Matlab has no option for Z=2048 in LTE toolbox)
    mac_meta_tx.oversampling = 2;                       % By how much do we oversample our ofdm packet compared to critical sampling (insert zeros at specturm edges before IFFT)?
    mac_meta_tx.codebook_index = 0;                     % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)
    mac_meta_tx.PLCF_type = 2;                          % Type 1 is 40 bits, Type 2 is 80 bits
    mac_meta_tx.rv = 0;                                 % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
    mac_meta_tx.network_id = de2bi(1e6,32,'left-msb');  % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY
    
    % temporary restrictions
    if mac_meta_tx.Z ~= 6144
        error('Z must be 6144.');
    end    
    
    % create tx
    verbose = 0;
    tx = dect_tx(verbose, mac_meta_tx);

    % additional rx configuration
    mac_meta_rx = mac_meta_tx;
    mac_meta_rx.N_RX = 1;

    % synchronization before the FFT (i.e. in time domain) based on the STF
    mac_meta_rx.synchronization.pre_FFT.active = false;
    if mac_meta_rx.synchronization.pre_FFT.active == true
    
        % symbol time offset (STO), i.e. detection, coarse peak search, fine peak search
        mac_meta_rx.synchronization.pre_FFT.sto_config = lib_rx.sync_STO_param(mac_meta_tx.u, mac_meta_tx.b, mac_meta_tx.oversampling);
        
        % carrier frequency offset (CFO), i.e. fractional and integer CFO
        mac_meta_rx.synchronization.pre_FFT.cfo_config = lib_rx.sync_CFO_param(mac_meta_tx.u);
        mac_meta_rx.synchronization.pre_FFT.cfo_config.active_fractional = false;
        mac_meta_rx.synchronization.pre_FFT.cfo_config.active_integer = false;
    end
    
    % synchronization in frequency domain based on STF and/or DRS
    mac_meta_rx.synchronization.post_FFT.sto_fractional = false;
    mac_meta_rx.synchronization.post_FFT.cfo_residual = false;

    % create rx
    rx = dect_rx(verbose, mac_meta_rx);

    % Local variables for a single MCS, required for parfor. Each worker can write into these arrays.
    snr_db_vec = snr_db_vec_global(cnt,:);
    
    % PCC
    n_bits_PCC_sent_local = zeros(1, numel(snr_db_vec));
    n_bits_PCC_error_local = zeros(1, numel(snr_db_vec));
    n_packets_PCC_sent_local = zeros(1, numel(snr_db_vec));
    n_packets_PCC_error_local = zeros(1, numel(snr_db_vec));
  
    % PDC
    n_bits_PDC_sent_local = zeros(1, numel(snr_db_vec));
    n_bits_PDC_error_local = zeros(1, numel(snr_db_vec));
    n_packets_PDC_sent_local = zeros(1, numel(snr_db_vec));
    n_packets_PDC_error_local = zeros(1, numel(snr_db_vec));

    %for i=1:numel(snr_db_vec)
    parfor i=1:numel(snr_db_vec)
        
        warning('off');
        
        % copy handle objects, changes within parfor are not permanent!
        txx = tx;
        rxx = rx;

        % run simulation over multiple packets
        result = simulate_packets(txx, rxx, snr_db_vec(i), n_packets_per_snr, max_harq_retransmissions);
        
        % each worker writes to local PCC result container
        n_bits_PCC_sent_local(1,i) = result.n_bits_PCC_sent;
        n_bits_PCC_error_local(1,i) = result.n_bits_PCC_error;
        n_packets_PCC_sent_local(1,i) = n_packets_per_snr;
        n_packets_PCC_error_local(1,i) = result.n_packets_PCC_error;         
        
        % each worker writes to local PDC result container
        n_bits_PDC_sent_local(1,i) = result.n_bits_PDC_sent;
        n_bits_PDC_error_local(1,i) = result.n_bits_PDC_error;
        n_packets_PDC_sent_local(1,i) = n_packets_per_snr;
        n_packets_PDC_error_local(1,i) = result.n_packets_PDC_error;
    end
    
    fprintf('Done! MCS %d of %d at %s\n', cnt, numel(mcs_index_vec), datestr(now,'HH:MM:SS'));
    
    % copy from local to global PCC results container
    n_bits_PCC_sent_global(cnt,:) = n_bits_PCC_sent_local;
    n_bits_PCC_error_global(cnt,:) = n_bits_PCC_error_local;
    n_packets_PCC_sent_global(cnt,:) = n_packets_PCC_sent_local;
    n_packets_PCC_error_global(cnt,:) = n_packets_PCC_error_local;    

    % copy from local to global PDC results container
    n_bits_PDC_sent_global(cnt,:) = n_bits_PDC_sent_local;
    n_bits_PDC_error_global(cnt,:) = n_bits_PDC_error_local;
    n_packets_PDC_sent_global(cnt,:) = n_packets_PDC_sent_local;
    n_packets_PDC_error_global(cnt,:) = n_packets_PDC_error_local;

    bps_global(cnt) = tx.phy_4_5.mcs.N_bps;
    tbs_global(cnt) = tx.phy_4_5.N_TB_bits;
    
    cnt = cnt + 1;
end

% save all variables
save('results/var_all.mat');

function [result] = simulate_packets(txx, rxx, snr_dB, n_packets_per_snr, max_harq_retransmissions)

    n_bits_PCC_sent = 0;
    n_bits_PCC_error = 0;
    n_packets_PCC_error = 0;

    n_bits_PDC_sent = 0;
    n_bits_PDC_error = 0;
    n_packets_PDC_error = 0;

    % how many antennas do we have?
    N_TX = txx.phy_4_5.tm_mode.N_TX;
    N_RX = rxx.mac_meta.N_RX;

    % create channel
    ch                      = lib_ch.rf_channel();
    ch.verbose              = 0;
    ch.verbose_cp           = txx.phy_4_5.numerology.N_b_CP*txx.mac_meta.oversampling;
    ch.type                 = 'rician';
    ch.amp                  = 1.0;
    ch.noise                = true;
    ch.snr_db               = snr_dB;
    ch.spectrum_occupied    = txx.phy_4_5.n_spectrum_occupied/txx.mac_meta.oversampling;
    ch.N_TX                 = N_TX;
    ch.N_RX                 = N_RX;
    ch.awgn_random_source   = 'global';
    ch.awgn_randomstream    = RandStream('mt19937ar','Seed', randi(1e9,[1 1]));
    ch.a_sto                = 0;
    ch.a_cfo                = 0;
    ch.a_err_phase          = 0;
    ch.r_random_source      = 'global';
    ch.r_seed               = randi(1e9,[1 1]);
    ch.r_sto                = 0;
    ch.r_cfo                = 0;
    ch.r_err_phase          = 0;
    ch.r_samp_rate          = txx.phy_4_5.numerology.B_u_b_DFT*txx.mac_meta.oversampling;
    ch.r_max_doppler        = 1.946;                            % 1.946 19.458
    ch.r_type               = 'TDL-v';
    ch.r_DS_desired         = 10^(-7.03 + 0.00*randn(1,1));
    ch.r_K                  = db2pow(9.0 + 0.00*randn(1,1));    %93e-9;
    ch.r_interpolation      = true;
    ch.r_gains_active       = true;
    ch.init_rayleigh_rician_channel();

    % adapt Wiener coefficients to channel conditions
    rxx.overwrite_wiener(1/10^(snr_dB/10), 20, 363e-9);

    % give rx handles so it can debug
    rxx.tx_handle = txx;
    rxx.ch_handle = ch;

    % how many bits does tx need?
    N_TB_bits = txx.phy_4_5.N_TB_bits;

    for j=1:1:n_packets_per_snr
        
        % generate random PCC bits
        if txx.mac_meta.PLCF_type == 1
            PCC_user_bits = randi([0 1], 40, 1);
        elseif txx.mac_meta.PLCF_type == 2
            PCC_user_bits = randi([0 1], 80, 1);
        end

        % generate bits
        PDC_user_bits = randi([0 1], N_TB_bits, 1);
        
        % harq abort conditions
        pcc_decoded_successfully = false;
        pdc_decoded_successfully = false;

        for z=0:1:max_harq_retransmissions

            % there is a specific order for the redundany version
            if mod(z,4) == 0
                txx.mac_meta.rv = 0; % initial transmission
                rxx.mac_meta.rv = 0; % initial transmission
            elseif mod(z,4) == 1
                txx.mac_meta.rv = 2;
                rxx.mac_meta.rv = 2;
            elseif mod(z,4) == 2
                txx.mac_meta.rv = 3;
                rxx.mac_meta.rv = 3;
            elseif mod(z,4) == 3
                txx.mac_meta.rv = 1;
                rxx.mac_meta.rv = 1;
            end

            % let tx create the packet
            samples_antenna_tx = txx.generate_packet(PCC_user_bits, PDC_user_bits);

            % pass samples through channel
            samples_antenna_rx = ch.pass_samples(samples_antenna_tx, 0);
            
            % make next channel impulse response independent from this one
            ch.reset_random_rayleigh_rician();

            % Now let rx decode the frame.
            % Rx can do so because it's mac_meta is the exact same.
            [PCC_user_bits_recovered, PDC_user_bits_recovered] = rxx.demod_decode_packet(samples_antenna_rx);
            
            % measure the BER uncoded
            
            n_bits_PCC_sent = n_bits_PCC_sent + numel(txx.packet_data.pcc_enc_dbg.d);
            n_bits_PCC_error = n_bits_PCC_error + sum(abs(double(txx.packet_data.pcc_enc_dbg.d) - double(rxx.packet_data.pcc_dec_dbg.d_hard)));                
            
            n_bits_PDC_sent = n_bits_PDC_sent + numel(txx.packet_data.pdc_enc_dbg.d);
            n_bits_PDC_error = n_bits_PDC_error + sum(abs(double(txx.packet_data.pdc_enc_dbg.d) - double(rxx.packet_data.pdc_dec_dbg.d_hard)));
            
            % we might be done
            if numel(PCC_user_bits_recovered) ~= 0
                pcc_decoded_successfully = true;
            end                

            % we might be done
            if numel(PDC_user_bits_recovered) ~= 0
                pdc_decoded_successfully = true;
            end
            
            % we continue sending retransmissions as long as not both were decoded correctly
            if pcc_decoded_successfully == true && pdc_decoded_successfully == true
                break;
            end
        end

        % delete harq buffer
        rxx.harq_buf_40 = [];
        rxx.harq_buf_80 = [];
        rxx.harq_buf = [];
        
        % check if frame was decoded correctly, maybe there's still an error despite all the harq iterations
        if pcc_decoded_successfully == false
            n_packets_PCC_error = n_packets_PCC_error + 1;
        end            

        % check if frame was decoded correctly, maybe there's still an error despite all the harq iterations
        if pdc_decoded_successfully == false
            n_packets_PDC_error = n_packets_PDC_error + 1;
        end
    end

    % delete channel
    delete(ch);

    result.n_bits_PCC_sent = n_bits_PCC_sent;
    result.n_bits_PCC_error = n_bits_PCC_error;
    result.n_packets_PCC_error = n_packets_PCC_error;

    result.n_bits_PDC_sent = n_bits_PDC_sent;
    result.n_bits_PDC_error = n_bits_PDC_error;
    result.n_packets_PDC_error = n_packets_PDC_error;
end