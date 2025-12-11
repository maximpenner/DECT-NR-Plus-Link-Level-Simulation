function [success] = single_packet_tx_rx()

    range.u_vec                 = [1,2];          % [1,2,4,8];
    range.b_vec                 = [1,2,4,8];      % [1,2,4,8,12,16];
    range.PacketLengthType_vec  = [0,1];
    range.PacketLength_vec      = 2:1:4;
    range.tm_mode_0_to_11_vec   = 0:1:0;
    range.mcs_index_vec         = 1:1:4;
    range.Z_vec                 = [2048,6144];
    range.oversampling_vec      = [1,2];
    range.codebook_index_vec    = 0:1:0;
    range.PLCF_type_vec         = [1,2];
    range.rv_vec                = 0:1:0;
    range.network_id_vec        = randi([1 2^32], 1, 1);
    range.verbosity             = 0;
    range.inner_repetitions     = 1;

    try
        dectnrp_regression.tx_config.parfor_harness(range, @test_per_config);
    catch ME
        s = dbstack;
        fprintf("Test %s, error message: %s\n\n", s(1).name, ME.message);
        success = false;
        return;
    end

    success = true;
end

function [] = test_per_config(tx_config)
    % create transmitter
    tx = dectnrp_tx.tx_t(tx_config);

    % create receiver, switch to equal weights to decrease simulation time
    rx_config = dectnrp_rx.rx_config_t();
    rx_config.channel_estimation_config.type = 'equal';
    rx = dectnrp_rx.rx_t(tx, rx_config);

    sync = dectnrp_sync.sync_t(tx);
    
    % create a DECT NR+ packet
    tx.generate_random_packet();
    
    % create channel configuration
    channel_config = dectnrp_channel.channel_config_t(tx_config.verbosity, 'AWGN', tx, rx);

    % change channel configuration
    channel_config.snr_db = 50;
    if tx_config.oversampling == 1
        channel_config.sto_fractional = 0;
    end

    % create channel
    channel = dectnrp_channel.channel_t(channel_config);
    
    % process one DECT NR+ packet
    samples_antenna_tx = tx.generate_random_packet();
    samples_antenna_ch = channel.pass_samples(samples_antenna_tx);
    samples_antenna_rx = sync.synchronize(samples_antenna_ch);
    rx.demod_decode_packet(samples_antenna_rx);
    
    % compare bits of transmitter and receiver
    assert(rx.are_tb_bits_equal(tx), "bits are not equal");
end