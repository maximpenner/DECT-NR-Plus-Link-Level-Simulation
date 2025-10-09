function [channel] = channel_example_factory(verbosity, type, tx, rx)

    % number of antennas at TX and RX
    N_TX = tx.phy_4_5.tm_mode.N_TX;
    N_RX = rx.rx_config.N_RX;

    % number of samples TX will generate
    n_samples_antenna_tx = tx.phy_4_5.n_packet_samples * tx.tx_config.oversampling;

    % RF channel parameters (see +lib_ch/channel.m) valid for all channel types.
    channel                     = lib_channel.channel_t();
    channel.verbosity           = verbosity;
    channel.type                = type;

    channel.N_TX                = N_TX;
    channel.N_RX                = N_RX;
    channel.spectrum_occupied   = tx.phy_4_5.n_spectrum_occupied/tx.tx_config.oversampling;

    channel.amp                 = 1.0;
    channel.sto_integer         = 123 + 2*n_samples_antenna_tx;
    channel.sto_fractional      = 0.36;
    channel.cfo                 = 1.7*(1/(tx.phy_4_5.numerology.N_b_DFT*tx.tx_config.oversampling));
    channel.err_phase           = deg2rad(123);
    channel.snr_db              = 30;
    
    if strcmp(channel.type, 'awgn')

        % nothing to do here

    elseif strcmp(channel.type, 'rayleigh') || strcmp(channel.type, 'rician')
    
        channel.r_samp_rate      = tx.phy_4_5.numerology.B_u_b_DFT*tx.tx_config.oversampling;
        channel.r_max_doppler    = 1.946;

        if strcmp(channel.type, 'rayleigh')
            channel.r_type       = 'TDL-iii';
        else
            channel.r_type       = 'TDL-iv';
        end

        channel.r_DS_desired     = 10^(-7.03);
        channel.r_K              = db2pow(9.0);

        channel.init_rayleigh_rician_channel();
    else
        assert(false, 'unknown channel type %d', channel.type);
    end
end
