clear all;
close all;

% This script illustrates how to use this DECT-2020 New Radio link-level simulation environment.
% It defines a transmitter and a receiver handle class object.
% Also, it initiates the wireless channel model (this step can be skipped and replaced by a custom channel model).
% The transmitter generates a packet, sends it through the wireless channel and finally the receiver decodes it.

rng('shuffle');

warning('on');

%% these variables need to be set before creating tx
mac_meta_tx.u = 1;                                  % mu = 1, 2, 4 or 8
mac_meta_tx.b = 1;                                  % beta = 1, 2, 4, 8, 12 or 16
mac_meta_tx.PacketLengthType = 0;                   % 0 for subslots, 1 for slots
mac_meta_tx.PacketLength = 4;                       % min is 1, max is 16 according to Table 6.2.1-2a in part 4
mac_meta_tx.tm_mode_0_to_11 = 0;                    % Table 7.2-1, mode determines whether transmission is closed loop or not, values range from 0 to 11
mac_meta_tx.mcs_index = 5;                          % Table A-1 in part 3, values range from 0 to 11
mac_meta_tx.Z = 6144;                               % 5.3, maximum codeblock size
mac_meta_tx.oversampling = 8;                       % By how much do we oversample our OFDM packet compared to critical sampling (insert zeros at spectrum edges before IFFT)?
mac_meta_tx.codebook_index = 0;                     % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)
mac_meta_tx.PLCF_type = 1;                          % Type 1 is 40 bits, Type 2 is 80 bits
mac_meta_tx.rv = 0;                                 % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
mac_meta_tx.network_id = de2bi(1e6,32,'left-msb');  % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY

% temporary restrictions
if mac_meta_tx.Z ~= 6144
    error('Z must be 6144.');
end

%% create tx
verbose = 2;                                        % show data during execution: 0 false, 1 only text, 2 text + plots
tx = dect_tx(verbose, mac_meta_tx);

%% generate tx signal

% generate random PCC bits
PCC_user_bits = [];
if mac_meta_tx.PLCF_type == 1
    PCC_user_bits = randi([0 1], 40, 1);
elseif mac_meta_tx.PLCF_type == 2
    PCC_user_bits = randi([0 1], 80, 1);
end

% how many PDC bits does tx need?
N_TB_bits = tx.phy_4_5.N_TB_bits;

% generate bits
PDC_user_bits = randi([0 1], N_TB_bits, 1);

% let tx create the packet
samples_antenna_tx = tx.generate_packet(PCC_user_bits, PDC_user_bits);

%% create rx

% assume the receiver has full knowledge of meta data at the transmitter (usually extracted from STF+PCC or blindly tested)
mac_meta_rx = mac_meta_tx;

% number of antennas at the receiver
mac_meta_rx.N_RX = 2;

% Synchronization before the FFT (i.e. in time domain) based on the STF:
%
% If synchronization before the FFT is turned on (i.e. mac_meta_rx.synchronization.pre_FFT.active = true), the receiver class
% dect_rx will try to synchronize a packet before decoding it. For that, the dect_rx method demod_decode_packet(samples_antenna_rx)
% must be called with samples_antenna_rx having more samples than samples_antenna_tx.
%
% If synchronization before the FFT is turned off (i.e. mac_meta_rx.synchronization.pre_FFT.active = false) the receiver class
% dect_rx will NOT try to synchronize a packet before decoding it. For that, the dect_rx method demod_decode_packet(samples_antenna_rx)
% must be called with samples_antenna_rx having the exact same number of samples as samples_antenna_tx.
%
mac_meta_rx.synchronization.pre_FFT.active = true;
if mac_meta_rx.synchronization.pre_FFT.active == true

    % symbol time offset (STO), i.e. detection, coarse peak search, fine peak search
    mac_meta_rx.synchronization.pre_FFT.sto_config = lib_rx.sync_STO_param(mac_meta_tx.u, mac_meta_tx.b, mac_meta_tx.oversampling);
    
    % carrier frequency offset (CFO), i.e. fractional and integer CFO
    mac_meta_rx.synchronization.pre_FFT.cfo_config = lib_rx.sync_CFO_param(mac_meta_tx.u);
    mac_meta_rx.synchronization.pre_FFT.cfo_config.active_fractional = true;
    mac_meta_rx.synchronization.pre_FFT.cfo_config.active_integer = true;
end

% synchronization in frequency domain based on STF and/or DRS
mac_meta_rx.synchronization.post_FFT.sto_fractional = true;
mac_meta_rx.synchronization.post_FFT.cfo_residual = true;

% create actual receiver
rx = dect_rx(verbose, mac_meta_rx);

%% create channel, can be replaced with a custom channel

% type can be awgn, rayleigh or rician
ch = lib_ch.rf_channel_example_factory('rayleigh', verbose, tx, rx, size(samples_antenna_tx, 1));

%% give rx handles so it can debug, e.g. perfect channel knowledge
rx.tx_handle = tx;
rx.ch_handle = ch;

%% pass tx signal through channel and yield the rx signal
samples_antenna_rx = ch.pass_samples(samples_antenna_tx, 0);

%% let rx decode the frame
[PCC_user_bits_recovered, PDC_user_bits_recovered] = rx.demod_decode_packet(samples_antenna_rx);
