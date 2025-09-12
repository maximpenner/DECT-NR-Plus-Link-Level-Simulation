function [ch] = rf_channel_example_factory(type, verbosity, tx, rx, n_samples_antenna_tx)

    % number of antennas at TX and RX
    N_TX = tx.phy_4_5.tm_mode.N_TX;
    N_RX = rx.config_rx.N_RX;

    % RF channel parameters (see +lib_rf_channel/rf_channel.m) valid for all channel types.
    ch                      = lib_ch.rf_channel();
    ch.verbosity            = verbosity;
    ch.verbosity_cp         = tx.phy_4_5.numerology.N_b_CP*tx.config_tx.oversampling;
    ch.type                 = type;
    ch.amp                  = 1.0;
    ch.noise                = true;
    ch.snr_db               = 30;
    ch.spectrum_occupied    = tx.phy_4_5.n_spectrum_occupied/tx.config_tx.oversampling;
    ch.N_TX                 = N_TX;
    ch.N_RX                 = N_RX;
    ch.awgn_random_source   = 'global';
    ch.awgn_randomstream    = RandStream('mt19937ar','Seed', randi(1e9,[1 1]));
    
    if strcmp(ch.type, 'awgn')
    
        % Parameters with a_ only used if ch.type = 'awgn'. It is AWGN with fixed values of STO, CFO and error phase.
        ch.a_sto            = 123 + 2*n_samples_antenna_tx;
        ch.a_cfo            = 1.7*(1/(tx.phy_4_5.numerology.N_b_DFT*tx.config_tx.oversampling));
        ch.a_err_phase      = deg2rad(123);
    
    elseif strcmp(ch.type, 'rayleigh') || strcmp(ch.type, 'rician')
    
        % Parameters with r_ only used if ch.type = 'rayleigh' or 'rician'.
        ch.r_random_source  = 'global';
        ch.r_seed           = randi(1e9,[1 1]);
        ch.r_sto            = 123 + 2*n_samples_antenna_tx;
        ch.r_cfo            = 1.7*(1/(tx.phy_4_5.numerology.N_b_DFT*tx.config_tx.oversampling));
        ch.r_err_phase      = deg2rad(123);
        ch.r_samp_rate      = tx.phy_4_5.numerology.B_u_b_DFT*tx.config_tx.oversampling;
        ch.r_max_doppler    = 1.946;                            % 1.946 19.458

        if strcmp(ch.type, 'rayleigh')
            ch.r_type       = 'TDL-iii';
        else
            ch.r_type       = 'TDL-iv';
        end

        ch.r_DS_desired     = 10^(-7.03 + 0.00*randn(1,1));
        ch.r_K              = db2pow(9.0 + 0.00*randn(1,1));    %93e-9;
        ch.r_interpolation  = true;
        ch.r_gains_active   = true;
        ch.init_rayleigh_rician_channel();
    else
        assert(false, 'unknown channel type %d', ch.type);
    end
end

