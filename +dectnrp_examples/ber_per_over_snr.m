function [] = ber_per_over_snr()
    % This script calculates BERs and PERs for DECT-2020 New Radio packets over various wireless channels (AWGN, Rayleigh, Rician) for different MCSs.
    % For each MCS and each SNR, the same number of packets is calculated.
    % Results are saved in the folder results/.
    % When this script is finished, the function ber_per_over_snr_plot.m can be used to plot the results.
    % Executing this script as is should take only a few seconds to minutes, depending of the system and multi-core capabilities (parfor).
    
    dectnrp_util.setup_script();
    
    if exist('results', 'dir')
        dectnrp_util.clear_directory('results');
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
    n_packets_per_snr = 100;
    
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
        tx_config = dectnrp_tx.tx_config_t();
    
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
        tx = dectnrp_tx.tx_t(tx_config);
    
        % create rx
        rx_config = dectnrp_rx.rx_config_t();
        rx_config.N_RX = 1;
        rx = dectnrp_rx.rx_t(tx, rx_config);
        
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
    
        %for i=1:numel(snr_db)
        parfor i=1:numel(snr_db)
            
            warning('off');
            
            % duplicate handle objects, changes within parfor change neither tx nor rx
            tx_duplicate = copy(tx);
            rx_duplate = copy(rx);
    
            % run simulation over multiple packets
            result = simulate_packets(tx_duplicate, rx_duplate, snr_db(i), n_packets_per_snr, harq_retransmissions);
            
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
    
        bps(cnt) = tx.tx_derived.mcs.N_bps;
        tbs(cnt) = tx.tx_derived.N_TB_bits;
    end
    
    % save all variables to file
    save('results/var_all.mat');
end

function [result] = simulate_packets(tx, rx, snr_dB, n_packets_per_snr, harq_retransmissions)

    n_bits_PCC_sent = 0;
    n_bits_PCC_error = 0;
    n_packets_PCC_error = 0;

    n_bits_PDC_sent = 0;
    n_bits_PDC_error = 0;
    n_packets_PDC_error = 0;

    % create channel configuration
    channel_config                      = dectnrp_channel.channel_config_t();
    channel_config.verbosity            = 0;
    channel_config.type                 = 'Rician';
    channel_config.N_TX                 = tx.tx_derived.tm_mode.N_TX;
    channel_config.N_RX                 = rx.rx_config.N_RX;
    channel_config.spectrum_occupied    = tx.tx_derived.n_spectrum_occupied/tx.tx_config.oversampling;
    channel_config.amp                  = 1.0;
    channel_config.sto_integer          = 0;
    channel_config.sto_fractional       = 0;
    channel_config.cfo                  = 0;
    channel_config.err_phase            = 0;
    channel_config.snr_db               = snr_dB;
    channel_config.r_samp_rate          = tx.tx_derived.numerology.B_u_b_DFT*tx.tx_config.oversampling;
    channel_config.r_max_doppler        = 1.946;
    channel_config.r_type               = 'TDL-v';
    channel_config.r_DS_desired         = 10^(-7.03 + 0.00*randn());
    channel_config.r_K                  = db2pow(9.0 + 0.00*randn());

    % create channel
    channel = dectnrp_channel.channel_t(channel_config);

    % adapt channel estimation weights to channel conditions
    rx.set_derived_weights(1/10^(snr_dB/10), 20, 363e-9);

    % how many bits does tx need?
    N_TB_bits = tx.tx_derived.N_TB_bits;

    for j=1:1:n_packets_per_snr
        
        % generate random PCC bits
        if tx.tx_config.PLCF_type == 1
            PCC_bits = randi([0 1], 40, 1);
        elseif tx.tx_config.PLCF_type == 2
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
                tx.tx_config.rv = 0;
                rx.tx_config.rv = 0;
            elseif mod(z,4) == 1
                tx.tx_config.rv = 2;
                rx.tx_config.rv = 2;
            elseif mod(z,4) == 2
                tx.tx_config.rv = 3;
                rx.tx_config.rv = 3;
            elseif mod(z,4) == 3
                tx.tx_config.rv = 1;
                rx.tx_config.rv = 1;
            end

            % let tx create the packet
            samples_antenna_tx = tx.generate_packet(PCC_bits, PDC_bits);

            % pass samples through channel
            samples_antenna_rx = channel.pass_samples(samples_antenna_tx, 0);
            
            % make next channel impulse response independent from this one
            channel.reset_random_Rayleigh_Rician();

            % now let rx decode the packet
            rx.demod_decode_packet(samples_antenna_rx);

            assert(numel(tx.packet_data.pcc_enc_dbg.d) == 196);
            
            % measure the BER uncoded
            n_bits_PCC_sent = n_bits_PCC_sent + numel(tx.packet_data.pcc_enc_dbg.d);
            n_bits_PCC_error = n_bits_PCC_error + rx.get_pcc_bit_errors_uncoded(tx);
            n_bits_PDC_sent = n_bits_PDC_sent + numel(tx.packet_data.pdc_enc_dbg.d);
            n_bits_PDC_error = n_bits_PDC_error + rx.get_pdc_bit_errors_uncoded(tx);
            
            % we might be done
            if rx.are_plcf_bits_equal(tx)
                pcc_decoded_successfully = true;
            end                

            % we might be done
            if rx.are_tb_bits_equal(tx)
                pdc_decoded_successfully = true;
            end
            
            % we continue sending re-transmissions as long as not both PCC and PDC are decoded correctly
            if pcc_decoded_successfully == true && pdc_decoded_successfully == true
                break;
            end
        end

        rx.clear_harq_buffers();
        
        % check if packet was decoded correctly, maybe there's still an error despite all the HARQ iterations
        if pcc_decoded_successfully == false
            n_packets_PCC_error = n_packets_PCC_error + 1;
        end            

        % check if packet was decoded correctly, maybe there's still an error despite all the HARQ iterations
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
