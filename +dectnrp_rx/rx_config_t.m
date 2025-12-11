classdef rx_config_t < matlab.mixin.Copyable
    
    properties
        % number of antennas at the receiver
        N_RX;

        % processing after the FFT in frequency domain
        sto_fractional_config;          % based on STF and optionally DRS
        sto_residual_config;            % based on DRS
        cfo_residual_config;            % based on DRS
        channel_estimation_config;      % based on optionally STF and DRS
        equalization_detection_config;
    end
    
    methods
        function obj = rx_config_t()
            obj = obj.set_example_values();
        end
    end

    methods (Hidden = true)
        function obj = set_example_values(obj)
            obj.N_RX = 2;

            % can be deactivated by leaving empty
            obj.sto_fractional_config = 1;
            obj.sto_residual_config = 1;
            obj.cfo_residual_config = 1;

            % Calculating Wiener filter coefficients takes very long for large packets.
            % For basic interpolation and for testing purposes, type can be set to 'equal' to speed up the process.
            obj.channel_estimation_config.type = 'Wiener';
            obj.channel_estimation_config.N_closest_DRS_pilots = 8;
            obj.channel_estimation_config.noise_estim = 1/10^(30/10);   % 30dB SNR
            obj.channel_estimation_config.f_d_hertz = 20;               % 20Hz Doppler
            obj.channel_estimation_config.tau_rms_sec = 363e-9;         % 363ns delay spread

            % can be deactivated by leaving empty
            obj.equalization_detection_config = 1;
        end
    end
end
