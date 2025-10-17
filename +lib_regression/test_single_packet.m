function [] = test_single_packet(config)
    % create transmitter
    tx = lib_tx.tx_t(config);
    rx = lib_rx.rx_t(tx);
    sync = lib_sync.sync_t(tx);
    
    % create a DECT NR+ packet
    tx.generate_random_packet();
    
    % create channel
    channel_config = lib_channel.config_t(config.verbosity, 'AWGN', tx, rx);
    channel_config.snr_db = 50;
    channel = lib_channel.channel_t(channel_config);
    
    % process one DECT NR+ packet
    samples_antenna_tx = tx.generate_random_packet();
    samples_antenna_ch = channel.pass_samples(samples_antenna_tx);
    samples_antenna_rx = sync.synchronize(samples_antenna_ch);
    rx.demod_decode_packet(samples_antenna_rx);
    
    % compare bits of transmitter and receiver
    assert(rx.are_tb_bits_equal(tx), "bits are not equal");
end