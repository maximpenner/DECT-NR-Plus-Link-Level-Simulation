classdef channel_config_t < matlab.mixin.Copyable
    
    properties
        verbosity;
        type;                   % AWGN, Rayleigh, Rician
        
        % variables used by all channel types
        N_TX;                   % number of transmit antennas
        N_RX;                   % number of receive antennas
        spectrum_occupied;      % share of spectrum occupied

        amp;                    % amplitude as linear factor, can be linked to large-scale pathloss
        sto_integer;            % symbol timing offset in samples, must be >= 0
        sto_fractional;         % symbol timing offset
        cfo;                    % carrier frequency offset in Hertz
        err_phase;              % error phase in radians
        snr_db;                 % in-band (i.e. occupied spectrum) signal-to-noise-ratio in dB

        % prefix r_ for Rayleigh and Rician
        r_samp_rate;            % system sample rate in Samples/s
        r_max_doppler;          % maximum doppler in Hertz
        r_type;                 % Exponential Decay: TDL-i, TDL-ii, TDL-iii(NLOS) or TDL-iv, TDL-v (LOS) etc.
        r_DS_desired;           % scaling factor of normalized delay spread (e.g. according to ITU-R M.2412-0, Table A1-43)
        r_K;                    % Rician fading K factor
    end
    
    methods
        function obj = channel_config_t(verbosity, type, tx, rx)

            if nargin > 0
                assert(~isempty(verbosity));
                assert(~isempty(type));
                assert(~isempty(tx));
                assert(~isempty(rx));

                obj = obj.set_example_values(verbosity, type, tx, rx);
            end
        end
    end

    methods (Hidden = true)
        function obj = set_example_values(obj, verbosity, type, tx, rx)
        
            % number of samples TX will generate
            n_samples_antenna_tx = tx.tx_derived.n_packet_samples * tx.tx_config.oversampling;
        
            % RF channel parameters (see +lib_ch/obj.m) valid for all channel types.
            obj.verbosity           = verbosity;
            obj.type                = type;
        
            obj.N_TX                = tx.tx_derived.tm_mode.N_TX;
            obj.N_RX                = rx.rx_config.N_RX;
            obj.spectrum_occupied   = tx.tx_derived.n_spectrum_occupied/tx.tx_config.oversampling;
        
            obj.amp                 = 1.0;
            obj.sto_integer         = 123 + 2*n_samples_antenna_tx;
            obj.sto_fractional      = 0.36;
            obj.cfo                 = 1.7*(1/(tx.tx_derived.numerology.N_b_DFT*tx.tx_config.oversampling));
            obj.err_phase           = deg2rad(123);
            obj.snr_db              = 30;
            
            if strcmp(obj.type, 'AWGN')
        
                % nothing to do here
        
            elseif strcmp(obj.type, 'Rayleigh') || strcmp(obj.type, 'Rician')
            
                obj.r_samp_rate      = tx.tx_derived.numerology.B_u_b_DFT*tx.tx_config.oversampling;
                obj.r_max_doppler    = 1.946;
        
                if strcmp(obj.type, 'Rayleigh')
                    obj.r_type       = 'TDL-iii';
                else
                    obj.r_type       = 'TDL-iv';
                end
        
                obj.r_DS_desired     = 10^(-7.03);
                obj.r_K              = db2pow(9.0);
            else
                assert(false, 'unknown channel type %d', obj.type);
            end
        end
    end
end
