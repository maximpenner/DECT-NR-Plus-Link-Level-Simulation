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
warning('off');

if exist('results', 'dir')
    lib_util.clear_directory('results');
else
    mkdir('results');
end

fprintf('Starting at %s\n', datestr(now,'HH:MM:SS'));

% choose mcs to simulate and maximum number of HARQ re-transmissions
mcs = [1,2,3,4];
harq_retransmissions = 0;

% simulation range for link level simulation
snr_db = -10:1.0:30;

% Packets per mcs and snr. Increase this number to get smoother curves.
n_packets_per_snr = 1.0e3;

% result container for PCC
n_bits_PCC_sent = zeros(numel(mcs), numel(snr_db(1,:)));        % BER uncoded
n_bits_PCC_error = zeros(numel(mcs), numel(snr_db(1,:)));       % BER uncoded
n_packets_PCC_sent = zeros(numel(mcs), numel(snr_db(1,:)));     % PER
n_packets_PCC_error = zeros(numel(mcs), numel(snr_db(1,:)));    % PER

% result container for PDC
n_bits_PDC_sent = zeros(numel(mcs), numel(snr_db(1,:)));        % BER uncoded
n_bits_PDC_error = zeros(numel(mcs), numel(snr_db(1,:)));       % BER uncoded
n_packets_PDC_sent = zeros(numel(mcs), numel(snr_db(1,:)));     % PER
n_packets_PDC_error = zeros(numel(mcs), numel(snr_db(1,:)));    % PER

% bits per symbol, transport block size
bps = zeros(numel(mcs), 1);
tbs = zeros(numel(mcs), 1);

for cnt = 1:numel(mcs)

    % configuration is initialized with exemplary values
    tx_config = lib_types.tx_config_t();

    % overwrite exemplary values
    tx_config.u = 1;
    tx_config.b = 1;
    tx_config.PacketLengthType = 0;
    tx_config.PacketLength = 2;
    tx_config.tm_mode_0_to_11 = 0;
    tx_config.mcs_index = mcs(cnt);
    tx_config.Z = 6144;
    tx_config.oversampling = 2;
    tx_config.codebook_index = 0;
    tx_config.PLCF_type = 2;
    tx_config.rv = 0;
    tx_config.network_id = de2bi(1e6,32,'left-msb');
    tx_config.verbosity = 0;
    
    % create tx
    tx = lib_types.tx_t(tx_config);

    % configuration is initialized with exemplary values
    rx_config = lib_types.rx_config_t(tx_config);

    % overwrite exemplary values
    rx_config.N_RX = 1;
    rx_config.pre_fft_config = [];

    % create rx
    rx = lib_types.rx_t(tx, rx_config);
    
    % PCC
    n_bits_PCC_sent_row = zeros(1, numel(snr_db));
    n_bits_PCC_error_row = zeros(1, numel(snr_db));
    n_packets_PCC_sent_row = zeros(1, numel(snr_db));
    n_packets_PCC_error_row = zeros(1, numel(snr_db));
  
    % PDC
    n_bits_PDC_sent_row = zeros(1, numel(snr_db));
    n_bits_PDC_error_row = zeros(1, numel(snr_db));
    n_packets_PDC_sent_row = zeros(1, numel(snr_db));
    n_packets_PDC_error_row = zeros(1, numel(snr_db));

    %for i=1:numel(snr_db_vec)
    parfor i=1:numel(snr_db)
        
        warning('off');
        
        % copy handle objects, changes within parfor are not permanent
        tx_cpy = copy(tx);
        rx_cpy = copy(rx);

        % run simulation over multiple packets
        result = simulate_packets(tx_cpy, rx_cpy, snr_db(i), n_packets_per_snr, harq_retransmissions);
        
        % each worker writes to local PCC result container
        n_bits_PCC_sent_row(1,i) = result.n_bits_PCC_sent;
        n_bits_PCC_error_row(1,i) = result.n_bits_PCC_error;
        n_packets_PCC_sent_row(1,i) = n_packets_per_snr;
        n_packets_PCC_error_row(1,i) = result.n_packets_PCC_error;         
        
        % each worker writes to local PDC result container
        n_bits_PDC_sent_row(1,i) = result.n_bits_PDC_sent;
        n_bits_PDC_error_row(1,i) = result.n_bits_PDC_error;
        n_packets_PDC_sent_row(1,i) = n_packets_per_snr;
        n_packets_PDC_error_row(1,i) = result.n_packets_PDC_error;
    end
    
    fprintf('Done! MCS %d of %d at %s\n', cnt, numel(mcs), datestr(now,'HH:MM:SS'));
    
    % copy from local to global PCC results container
    n_bits_PCC_sent(cnt,:) = n_bits_PCC_sent_row;
    n_bits_PCC_error(cnt,:) = n_bits_PCC_error_row;
    n_packets_PCC_sent(cnt,:) = n_packets_PCC_sent_row;
    n_packets_PCC_error(cnt,:) = n_packets_PCC_error_row;    

    % copy from local to global PDC results container
    n_bits_PDC_sent(cnt,:) = n_bits_PDC_sent_row;
    n_bits_PDC_error(cnt,:) = n_bits_PDC_error_row;
    n_packets_PDC_sent(cnt,:) = n_packets_PDC_sent_row;
    n_packets_PDC_error(cnt,:) = n_packets_PDC_error_row;

    bps(cnt) = tx.phy_4_5.mcs.N_bps;
    tbs(cnt) = tx.phy_4_5.N_TB_bits;
end

% save all variables to file
save('results/var_all.mat');

function [result] = simulate_packets(tx_cpy, rx_cpy, snr_dB, n_packets_per_snr, harq_retransmissions)

    n_bits_PCC_sent = 0;
    n_bits_PCC_error = 0;
    n_packets_PCC_error = 0;

    n_bits_PDC_sent = 0;
    n_bits_PDC_error = 0;
    n_packets_PDC_error = 0;

    % how many antennas do we have?
    N_TX = tx_cpy.phy_4_5.tm_mode.N_TX;
    N_RX = rx_cpy.rx_config.N_RX;

    % create channel
    channel                     = lib_channel.channel_t();
    channel.verbosity           = 0;
    channel.type                = 'rician';
    channel.N_TX                = N_TX;
    channel.N_RX                = N_RX;
    channel.spectrum_occupied   = tx_cpy.phy_4_5.n_spectrum_occupied/tx_cpy.tx_config.oversampling;
    channel.amp                 = 1.0;
    channel.sto_integer         = 0;
    channel.sto_fractional      = 0;
    channel.cfo                 = 0;
    channel.err_phase           = 0;
    channel.snr_db              = snr_dB;
    channel.r_samp_rate         = tx_cpy.phy_4_5.numerology.B_u_b_DFT*tx_cpy.tx_config.oversampling;
    channel.r_max_doppler       = 1.946;
    channel.r_type              = 'TDL-v';
    channel.r_DS_desired        = 10^(-7.03 + 0.00*randn(1,1));
    channel.r_K                 = db2pow(9.0 + 0.00*randn(1,1));
    channel.init_rayleigh_rician_channel();

    % adapt Wiener coefficients to channel conditions
    rx_cpy.set_wiener(1/10^(snr_dB/10), 20, 363e-9);

    % how many bits does tx need?
    N_TB_bits = tx_cpy.phy_4_5.N_TB_bits;

    for j=1:1:n_packets_per_snr
        
        % generate random PCC bits
        if tx_cpy.tx_config.PLCF_type == 1
            PCC_bits = randi([0 1], 40, 1);
        elseif tx_cpy.tx_config.PLCF_type == 2
            PCC_bits = randi([0 1], 80, 1);
        end

        % generate bits
        PDC_bits = randi([0 1], N_TB_bits, 1);
        
        % HARQ abort conditions
        pcc_decoded_successfully = false;
        pdc_decoded_successfully = false;

        for z=0:1:harq_retransmissions

            % there is a specific order for the redundancy version
            if mod(z,4) == 0
                tx_cpy.tx_config.rv = 0;
                rx_cpy.tx_config.rv = 0;
            elseif mod(z,4) == 1
                tx_cpy.tx_config.rv = 2;
                rx_cpy.tx_config.rv = 2;
            elseif mod(z,4) == 2
                tx_cpy.tx_config.rv = 3;
                rx_cpy.tx_config.rv = 3;
            elseif mod(z,4) == 3
                tx_cpy.tx_config.rv = 1;
                rx_cpy.tx_config.rv = 1;
            end

            % let tx create the packet
            samples_antenna_tx = tx_cpy.generate_packet(PCC_bits, PDC_bits);

            % pass samples through channel
            samples_antenna_rx = channel.pass_samples(samples_antenna_tx, 0);
            
            % make next channel impulse response independent from this one
            channel.reset_random_rayleigh_rician();

            % now let rx decode the frame
            rx_cpy.demod_decode_packet(samples_antenna_rx);

            assert(numel(tx_cpy.packet_data.pcc_enc_dbg.d) == 196);
            
            % measure the BER uncoded
            n_bits_PCC_sent = n_bits_PCC_sent + numel(tx_cpy.packet_data.pcc_enc_dbg.d);
            n_bits_PCC_error = n_bits_PCC_error + rx_cpy.get_pcc_bit_errors_uncoded(tx_cpy);
            n_bits_PDC_sent = n_bits_PDC_sent + numel(tx_cpy.packet_data.pdc_enc_dbg.d);
            n_bits_PDC_error = n_bits_PDC_error + rx_cpy.get_pdc_bit_errors_uncoded(tx_cpy);
            
            % we might be done
            if rx_cpy.are_plcf_bits_equal(tx_cpy)
                pcc_decoded_successfully = true;
            end                

            % we might be done
            if rx_cpy.are_tb_bits_equal(tx_cpy)
                pdc_decoded_successfully = true;
            end
            
            % we continue sending re-transmissions as long as not both PCC and PDC are decoded correctly
            if pcc_decoded_successfully == true && pdc_decoded_successfully == true
                break;
            end
        end

        rx_cpy.clear_harq_buffers();
        
        % check if frame was decoded correctly, maybe there's still an error despite all the HARQ iterations
        if pcc_decoded_successfully == false
            n_packets_PCC_error = n_packets_PCC_error + 1;
        end            

        % check if frame was decoded correctly, maybe there's still an error despite all the HARQ iterations
        if pdc_decoded_successfully == false
            n_packets_PDC_error = n_packets_PDC_error + 1;
        end
    end

    delete(channel);

    result.n_bits_PCC_sent = n_bits_PCC_sent;
    result.n_bits_PCC_error = n_bits_PCC_error;
    result.n_packets_PCC_error = n_packets_PCC_error;

    result.n_bits_PDC_sent = n_bits_PDC_sent;
    result.n_bits_PDC_error = n_bits_PDC_error;
    result.n_packets_PDC_error = n_packets_PDC_error;
end
