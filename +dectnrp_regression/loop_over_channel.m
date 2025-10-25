function [] = loop_over_channel()

    for i=1:1:100

        % create channel configuration
        channel_config                      = dectnrp_channel.config_t();
        channel_config.verbosity            = 0;
        channel_config.type                 = 'Rayleigh';
        channel_config.N_TX                 = randi([1 8]);
        channel_config.N_RX                 = randi([1 8]);
        channel_config.spectrum_occupied    = 0.125 + rand()*(0.99-0.125);
        channel_config.amp                  = 0.125 + rand()*(0.99-0.125);
        channel_config.sto_integer          = randi([0 10e3]);
        channel_config.sto_fractional       = 0;
        channel_config.cfo                  = 0;
        channel_config.err_phase            = 0;
        channel_config.snr_db               = 1000;
        channel_config.r_samp_rate          = randi([1e6 10e6]);
        channel_config.r_max_doppler        = 0;
        channel_config.r_type               = 'TDL-iii';
        channel_config.r_DS_desired         = 10^(-7.03 + 0.50*randn());
        channel_config.r_K                  = db2pow(9.0 + 0.50*randn());
    
        % create channel
        channel = dectnrp_channel.channel_t(channel_config);
    
        % create random noise signal
        samples_antenna_tx = randn(10e3, channel_config.N_TX);

        % create a random TX time
        channel_time_in_seconds = randi([0 1e6]) + rand();
    
        % generate the same samples at the same time
        samples_antenna_ch_A = channel.pass_samples(samples_antenna_tx, channel_time_in_seconds);
        channel.reset_random_Rayleigh_Rician();
        samples_antenna_ch_B = channel.pass_samples(samples_antenna_tx, channel_time_in_seconds);

        assert(max(abs(samples_antenna_ch_A - samples_antenna_ch_B), [], 'all') < 1e-6);
    end
end