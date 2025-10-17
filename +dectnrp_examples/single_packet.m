function [] = single_packet()
    clear all;
    close all;
    
    % This script illustrates how to use this DECT-2020 New Radio link-level simulation environment.
    
    rng('shuffle');
    warning('on');

    % The type config_t contains the minimal set of variables that define the structure and size of a DECT NR+ packet.
    % It is initialized with sample values that can be overwritten.
    config = dectnrp_tx.config_t();
    
    % create transmitter, receiver and synchronization
    tx = dectnrp_tx.tx_t(config);
    rx = dectnrp_rx.rx_t(tx);
    sync = dectnrp_sync.sync_t(tx);
    
    % create channel
    channel_config = dectnrp_channel.config_t(config.verbosity, 'Rayleigh', tx, rx);
    channel = dectnrp_channel.channel_t(channel_config);
    
    % process one DECT NR+ packet
    samples_antenna_tx = tx.generate_random_packet();
    samples_antenna_ch = channel.pass_samples(samples_antenna_tx);
    samples_antenna_rx = sync.synchronize(samples_antenna_ch);
    rx.demod_decode_packet(samples_antenna_rx);
    
    % compare bits of transmitter and receiver
    if rx.are_tb_bits_equal(tx)
        fprintf('Packet decoded correctly with an SINR of %f dB\n', rx.get_sinr(tx));
    else
        fprintf('Packet decoded incorrectly\n');
    end
end
