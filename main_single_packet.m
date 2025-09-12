clear all;
close all;

% This script illustrates how to use this DECT-2020 New Radio link-level simulation environment.

rng('shuffle');
warning('on');

% configurations are initialized with exemplary values
tx_config = tx_config_t();
rx_config = rx_config_t(tx_config);

% create transmitter, receiver and channel
tx = tx_t(tx_config);
rx = rx_t(tx, rx_config);
channel = lib_channel.channel_example_factory(tx_config.verbosity, 'rayleigh', tx, rx);

% create a DECT NR+ packet
samples_antenna_tx = tx.generate_random_packet();

% send it through the wireless channel
samples_antenna_rx = channel.pass_samples(samples_antenna_tx);

% let the receiver demodulate and decode it
[plcf_bits_recovered, tb_bits_recovered] = rx.demod_decode_packet(samples_antenna_rx);

% show result
if rx.are_tb_bits_equal(tx)
    fprintf('Packet decoded correctly with an SINR of %f dB\n', rx.get_sinr(tx));
else
    fprintf('Packet decoded incorrectly\n');
end
