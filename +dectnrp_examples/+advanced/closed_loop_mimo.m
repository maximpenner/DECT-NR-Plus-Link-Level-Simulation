function [] = closed_loop_mimo()


    dectnrp_util.setup_script();

    % In order to transmit packets via closed loop mimo, feedback by the
    % receiver is required. In this case the receiver does an estimation on
    % how many spatial streams and which precoding matrix works best.
    % In 5G there are specific reference symbols (CSI-RS), which are added 
    % after precoding, so the receiver can always estimate the true H 
    % regardless if precoding is used or not. In DECT-2020 NR we dont have
    % these reference symbols, instead there is a "channel sounding packet"
    % defined in TS 103 636-3 p. 31, which states that the "Channel
    % sounding packet shall always be transmitted with identity precoding
    % matrix (and open loop mimo)". So in order to activate closed loop
    % mimo transmissions, regular "Channel sounding packets" have to be
    % scheduled which are send with codebook index 0 and open loop MIMO. 
    % (So the receiver will estimate H and not W*H).

    % first packet is send with W = I (channel sounding packet)
    tx_config = dectnrp_tx.tx_config_t();
    tx_config.b = 2;
    tx_config.mcs_index = 1;
    tx_config.verbosity = 0;
    tx_config.PacketLengthType = 1;
    tx_config.PacketLength = 2;
    tx_config.tm_mode_0_to_11 = 9;
    tx_config.oversampling = 1;

    % open loop MIMO requires N_RX = N_TX
    rx_config = dectnrp_rx.rx_config_t();
    switch tx_config.tm_mode_0_to_11
        case {0}
            rx_config.N_RX = 1;
        case {1,2,3,4}
            rx_config.N_RX = 2;
        case {5,6,7,8,9}
            rx_config.N_RX = 4;
        case {10,11}
            rx_config.N_RX = 8;
    end
    
    % create transmitter, receiver and synchronization
    tx = dectnrp_tx.tx_t(tx_config);
    rx = dectnrp_rx.rx_t(tx, rx_config);

    % SNR dB vec
    snr_db_vec = 2:2:30;

    % packets per snr
    n_packets = 50;

    % result container
    sinr_1 = zeros(1, numel(snr_db_vec));
    sinr_2 = zeros(1, numel(snr_db_vec));
    ri = zeros(1, numel(snr_db_vec));
    
    % create channel configuration
    channel_config                      = dectnrp_channel.channel_config_t();
    channel_config.verbosity            = 0;
    channel_config.type                 = 'Rayleigh';
    channel_config.N_TX                 = tx.tx_derived.tm_mode.N_TX;
    channel_config.N_RX                 = rx.rx_config.N_RX;
    channel_config.spectrum_occupied    = tx.tx_derived.n_spectrum_occupied/tx.tx_config.oversampling;
    channel_config.amp                  = 1.0;
    channel_config.sto_integer          = 0;
    channel_config.sto_fractional       = 0;
    channel_config.cfo                  = 0;
    channel_config.err_phase            = 0;
    channel_config.snr_db               = snr_db_vec(end);
    channel_config.r_samp_rate          = tx.tx_derived.numerology.B_u_b_DFT*tx.tx_config.oversampling;
    channel_config.r_max_doppler        = 1.946;
    channel_config.r_type               = 'TDL-iii';
    channel_config.r_DS_desired         = 0;
    channel_config.r_K                  = db2pow(9.0);

    % create new channel
    channel = dectnrp_channel.channel_t(channel_config);
    channel.set_randomstream_for_channel_reproducibility(randi([1 1e9]));

    parfor i=1:numel(snr_db_vec)
        [sinr_1(i), sinr_2(i), ri(i)] = run_closed_loop(n_packets, snr_db_vec(i), tx, rx, channel);
    end

    figure;
    tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');
    
    nexttile;
    plot(snr_db_vec, sinr_1, 'LineWidth', 1.5);
    grid on;
    ylabel('SINR (dB)');
    title('Channel Sounding SINR');
    
    nexttile;
    plot(snr_db_vec, sinr_2, 'LineWidth', 1.5);
    grid on;
    ylabel('SINR (dB)');
    title('Closed-Loop SINR');
    
    nexttile;
    plot(snr_db_vec, ri, 'LineWidth', 1.5);
    grid on;
    xlabel('SNR (dB)');
    ylabel('RI');
    title('Average Rank Indicator');
    ylim([0 5])
    

end

function [sinr_channel_sounding_avg, sinr_cl_avg, ri_avg] = run_closed_loop(n_packets, snr, tx, rx, channel)
    
    % overwrite snr in channel
    channel.channel_config.snr_db = snr;

    % create containers
    sinr_channel_sounding_avg = 0;
    sinr_cl_avg = 0;
    ri_avg = 0;

    % create tx and rx objects for closed loop
    tx_cl_config = tx.tx_config.copy;
    rx_cl_config = rx.rx_config.copy;

    % loop over n_packets
    for j=1:n_packets

        % copy objects
        txx = tx.copy;
        rxx = rx.copy;

        % process one DECT NR+ packet
        samples_antenna_tx = txx.generate_random_packet();
        samples_antenna_ch = channel.pass_samples(samples_antenna_tx);
        rxx.demod_decode_packet(samples_antenna_ch);
    
        if rxx.are_tb_bits_equal(txx)
                sinr_channel_sounding = rxx.get_sinr(txx);
        else
                sinr_channel_sounding = -100;
        end

        ri = rxx.rx_derived.mimo.N_TS;
        
        % translate the receiver feedback to a transmission mode
        switch txx.tx_config.tm_mode_0_to_11
            case {2, 4}
                switch rxx.rx_derived.mimo.N_TS
                    case 1
                        closed_loop_tm_mode_0_to_11 = 3;
                    case 2
                        closed_loop_tm_mode_0_to_11 = 4;
                    otherwise
                        error("something went wrong")
                end
            case {6, 9}
                 switch rxx.rx_derived.mimo.N_TS
                    case 1
                        closed_loop_tm_mode_0_to_11 = 7;
                    case 2
                        closed_loop_tm_mode_0_to_11 = 8;
                    case 4
                        closed_loop_tm_mode_0_to_11 = 9;                    
                    otherwise
                        error("something went wrong")
                end
            case 11
                closed_loop_tm_mode_0_to_11 = 11;
                warning("only open loop mimo possible for N_TX = 8")
            otherwise
                error("something went wrong")
        
        end
    
        % create new transmitter and feed in feedback
        tx_cl_config.tm_mode_0_to_11 = closed_loop_tm_mode_0_to_11;
        tx_cl_config.codebook_index = rxx.rx_derived.mimo.codebook_index;
        txx_cl = dectnrp_tx.tx_t(tx_cl_config);
    
        %rxx_config = rxx.rx_config;
        rxx_cl = dectnrp_rx.rx_t(txx_cl, rx_cl_config);
    
        samples_antenna_txx = txx_cl.generate_random_packet();
        samples_antenna_chh = channel.pass_samples(samples_antenna_txx, 5 * 10^(-3));
        %channel.reset_random_Rayleigh_Rician();
        %samples_antenna_chh = channel.pass_samples(samples_antenna_txx);
        %samples_antenna_rxx = sync_new.synchronize(samples_antenna_chh);
        rxx_cl.demod_decode_packet(samples_antenna_chh);
    
        if rxx_cl.are_tb_bits_equal(txx_cl)
                sinr_cl = rxx_cl.get_sinr(txx_cl);
        else
                sinr_cl = -100;
        end

        channel.reset_random_Rayleigh_Rician();

        % save results
        sinr_channel_sounding_avg = sinr_channel_sounding_avg + sinr_channel_sounding;
        sinr_cl_avg = sinr_cl_avg + sinr_cl;
        ri_avg = ri_avg + ri;
    end

    % average results
    sinr_channel_sounding_avg = sinr_channel_sounding_avg/n_packets;
    sinr_cl_avg = sinr_cl_avg/n_packets;
    ri_avg = ri_avg/n_packets;
end
