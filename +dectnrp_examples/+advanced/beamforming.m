function [] = beamforming()

    % create packet configuration with four TX antennas, single spatial and single transmit stream
    config = dectnrp_tx.config_t();
    config.tm_mode_0_to_11 = 7;
    config.mcs_index = 1;
    config.verbosity = 0;

    % 27 different beamforming matrices for mode 7
    codebook_index_max = 27;
    codebook_index_vec = 0:codebook_index_max;
    
    % create transmitter
    tx = dectnrp_tx.tx_t(config);

    % create receiver
    rx = dectnrp_rx.rx_t(tx);
    rx.N_RX = 4;

    % create channel configuration
    channel_config                      = dectnrp_channel.config_t();
    channel_config.verbosity            = 0;
    channel_config.type                 = 'Rayleigh';
    channel_config.N_TX                 = tx.derived.tm_mode.N_TX;
    channel_config.N_RX                 = rx.N_RX;
    channel_config.spectrum_occupied    = tx.derived.n_spectrum_occupied/tx.config.oversampling;
    channel_config.amp                  = 1.0;
    channel_config.sto_integer          = 0;
    channel_config.sto_fractional       = 0;
    channel_config.cfo                  = 0;
    channel_config.err_phase            = 0;
    channel_config.snr_db               = 20;
    channel_config.r_samp_rate          = tx.derived.numerology.B_u_b_DFT*tx.config.oversampling;
    channel_config.r_max_doppler        = 1.946;
    channel_config.r_type               = 'TDL-iii';
    channel_config.r_DS_desired         = 0;
    channel_config.r_K                  = db2pow(9.0);

    % create new channel
    channel = dectnrp_channel.channel_t(channel_config);
    channel.set_randomstream_for_channel_reproducibility(randi([1 1e9]));

    % number of packets we evaluate
    N_packets = 100;

    % final results container
    sinr_improvement = zeros(N_packets, 1);

    % number of packets per beamforming matrix
    for n = 1:1:N_packets
        sinr_improvement(n) = run_single_codebook_index(tx, channel, rx, codebook_index_vec);
    end

    % plot results
    figure(1)
    clf()
    hold on
    plot(sinr_improvement);
    yline(pow2db(mean(db2pow(sinr_improvement))))
    grid on
    legend(["SINR improvement in dB relative to non-beamformed packet with optimal codebook index", "average across all packets"]);
    ylim([-20 20])
    ylabel("SINR Improvement in dB")
    xlabel("Packet Index");
end

function [sinr] = run_single_codebook_index(tx, channel, rx, codebook_index_vec)

    sinr = zeros(numel(codebook_index_vec), 1);

    % create a packet without beamforming
    tx.config.codebook_index = 0;
    samples_antenna_tx = tx.generate_random_packet();

    % simulate packet without beamforming until it was decoded successfully
    while sinr(1) <= 0
        % determine some random channel time
        channel_time = randi([0 1e6]) + rand();
    
        % pass packet without beamforming through channel
        channel.reset_random_Rayleigh_Rician();
        samples_antenna_ch = channel.pass_samples(samples_antenna_tx, channel_time);
        
        % let RX decode the packet without beamforming
        rx.clear_harq_buffers();
        rx.demod_decode_packet(samples_antenna_ch);
    
        % measurement without beamforming
        if rx.are_tb_bits_equal(tx)
            sinr(1) = rx.get_sinr(tx);
        else
            sinr(1) = -100;
        end
    end

    % start with the second index
    for i = 2:numel(codebook_index_vec)
        
        % create same packet with beamforming
        tx.config.codebook_index = codebook_index_vec(i);
        samples_antenna_tx = tx.generate_packet(tx.packet_data.plcf_bits, tx.packet_data.tb_bits);
        
        % pass beamformed signal through the exact same channel at the exact same time
        channel.reset_random_Rayleigh_Rician();
        samples_antenna_ch = channel.pass_samples(samples_antenna_tx, channel_time);

        % let RX decode the packet with beamforming
        rx.clear_harq_buffers();
        rx.demod_decode_packet(samples_antenna_ch);

        % measurement with beamforming
        if rx.are_tb_bits_equal(tx)
            sinr(i) = rx.get_sinr(tx);
        else
            sinr(i) = -100;
        end
    end

    % SINR improvement relative to packet without beamforming
    sinr = sinr(2:end) - sinr(1);
    sinr = max(sinr);
end