classdef rx_config_t < matlab.mixin.Copyable
    
    properties
        % number of antennas at the receiver
        N_RX

        % structure describes the signal processing based on the STF prior to the FFT
        pre_fft_config

        % synchronization after the FFT in frequency domain
        post_fft_sto_fractional_config  % based on STF and optionally DRS
        post_FFT_sto_residual_config    % based on DRS
        post_FFT_cfo_residual_config    % based on DRS

        channel_estimation_config

        equalization_detection_config

        % time and frequency domain templates
        stf_templates
    end
    
    methods
        function obj = rx_config_t(tx_config)
            obj = obj.set_example_values(tx_config);
            assert(obj.is_valid());
        end

        function ret = is_valid(obj)
            ret = true;

            if isempty(obj.stf_templates)
                ret = false;
            end
        end
    end

    methods (Hidden = true)
        function obj = set_example_values(obj, tx_config)
            obj.N_RX = 2;

            obj.pre_fft_config = lib_rx.lib_pre_fft.config(tx_config.u, tx_config.b, tx_config.oversampling);

            obj.stf_templates = lib_rx.stf_templates(tx_config);

            % can be deactivated by leaving empty
            obj.post_fft_sto_fractional_config = 1;
            obj.post_FFT_sto_residual_config = 1;
            obj.post_FFT_cfo_residual_config = 1;

            obj.channel_estimation_config.N_closest_DRS_pilots = 8;
            obj.channel_estimation_config.noise_estim = 1/10^(30/10);   % 30dB SNR
            obj.channel_estimation_config.f_d_hertz = 20;               % 20Hz Doppler
            obj.channel_estimation_config.tau_rms_sec = 363e-9;         % 363ns delay spread

            % can be deactivated by leaving empty
            obj.equalization_detection_config = 1;
        end
    end
end
