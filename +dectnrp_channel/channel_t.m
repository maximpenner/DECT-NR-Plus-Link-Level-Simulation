classdef channel_t < handle
    
    properties
        config;
        
        % reference to the matlab object
        r_matlab_MIMO_obj;

        % appended samples by Rayleigh channel
        r_appendix;
    end
    
    methods (Static = true, Access = public)
       
       % Init parameters outside of constructor, or use the rf_channel_example_factory.
       % This is typically more convenient than controlling all parameters through class methods.
       function obj = channel_t(config)
            assert(isa(config, "dectnrp_channel.config_t"));
            assert(config.is_valid());

            obj.config = config;

            if strcmp(obj.config.type,'Rayleigh') || strcmp(obj.config.type,'Rician')
                obj.init_Rayleigh_Rician_channel();
            end
       end
    end
    
    methods
        function init_Rayleigh_Rician_channel(obj)
            
            % lookup PDP
            [pathDelays_beforeInterpolation, avgPathGains_beforeInterpolation] = dectnrp_channel.get_PDP_from_literature(obj.config.r_type, ...
                                                                                                                         obj.config.r_DS_desired);
            
            % interpolate to system sample rate
            [pathDelays, avgPathGains] = obj.interpolate_power_delay_profile_to_system_sample_rate(pathDelays_beforeInterpolation, ...
                                                                                                   avgPathGains_beforeInterpolation);
            
            assert(numel(pathDelays) == numel(unique(pathDelays)));
            assert(sum(pathDelays < 0) == 0);
            assert(sum(isnan(avgPathGains)) == 0);
            assert(sum(~isfinite(avgPathGains)) == 0);
            
            % determine how many samples the channel is adding
            Ts = 1/obj.config.r_samp_rate;
            obj.r_appendix = ceil(max(pathDelays)/Ts);

            % create Matlab channel object
            if strcmp(obj.config.type, 'Rayleigh') == true
                
                assert(strcmp(obj.config.r_type, 'TDL-i') || ...
                       strcmp(obj.config.r_type, 'TDL-ii') || ...
                       strcmp(obj.config.r_type, 'TDL-iii') || ...
                       strcmp(obj.config.r_type, 'InH NLOS') || ...
                       strcmp(obj.config.r_type, 'UMi NLOS') || ...
                       strcmp(obj.config.r_type, 'UMi OtI'));

                obj.r_matlab_MIMO_obj = comm.MIMOChannel('SampleRate', obj.config.r_samp_rate, ...
                                                         'PathDelays', pathDelays, ...
                                                         'AveragePathGains', avgPathGains, ...
                                                         'NormalizePathGains', true, ...
                                                         'FadingDistribution', 'Rayleigh', ... % KFactor, DirectPathDopplerShift, DirectPathInitialPhase
                                                         'MaximumDopplerShift', obj.config.r_max_doppler, ...
                                                         'DopplerSpectrum', doppler('Jakes'), ...
                                                         'SpatialCorrelationSpecification', 'Separate Tx Rx', ... %'NumTransmitAntennas', obj.config.N_TX, 'NumReceiveAntennas', obj.config.N_RX, ...
                                                         'TransmitCorrelationMatrix', eye(obj.config.N_TX), ...
                                                         'ReceiveCorrelationMatrix', eye(obj.config.N_RX), ... % SpatialCorrelationMatrix
                                                         'AntennaSelection', 'off', ...
                                                         'NormalizeChannelOutputs', false, ...
                                                         'FadingTechnique', 'Sum of sinusoids', ...
                                                         'NumSinusoids', 48, ...
                                                         'InitialTimeSource', 'Input Port', ... % InitialTime
                                                         'RandomStream', 'Global stream', ... % Seed
                                                         'PathGainsOutputPort', false, ...
                                                         'Visualization', 'off'); %AntennaPairsToDisplay, PathsForDopplerDisplay, SamplesToDisplay
                
            elseif strcmp(obj.config.type, 'Rician') == true
                
                assert(strcmp(obj.config.r_type, 'TDL-iv') || ...
                       strcmp(obj.config.r_type, 'TDL-v') || ...
                       strcmp(obj.config.r_type, 'InH LOS') || ...
                       strcmp(obj.config.r_type, 'UMi LOS'));
                
                obj.r_matlab_MIMO_obj = comm.MIMOChannel('SampleRate', obj.config.r_samp_rate, ...
                                                         'PathDelays', pathDelays, ...
                                                         'AveragePathGains', avgPathGains, ...
                                                         'NormalizePathGains', true, ...
                                                         'FadingDistribution', 'Rician', ...
                                                         'KFactor', obj.config.r_K, ...
                                                         'DirectPathDopplerShift',obj.config.r_max_doppler*0.7, ...
                                                         'DirectPathInitialPhase',0, ...
                                                         'MaximumDopplerShift', obj.config.r_max_doppler, ...
                                                         'DopplerSpectrum', doppler('Jakes'), ...
                                                         'SpatialCorrelationSpecification', 'Separate Tx Rx', ... %'NumTransmitAntennas', obj.config.N_TX, 'NumReceiveAntennas', obj.config.N_RX, ...
                                                         'TransmitCorrelationMatrix', eye(obj.config.N_TX), ...
                                                         'ReceiveCorrelationMatrix', eye(obj.config.N_RX), ... % SpatialCorrelationMatrix
                                                         'AntennaSelection', 'off', ...
                                                         'NormalizeChannelOutputs', false, ...
                                                         'FadingTechnique', 'Sum of sinusoids', ...
                                                         'NumSinusoids', 48, ...
                                                         'InitialTimeSource', 'Input Port', ... % InitialTime
                                                         'RandomStream', 'Global stream', ... % Seed
                                                         'PathGainsOutputPort', false, ...
                                                         'Visualization', 'off'); %AntennaPairsToDisplay, PathsForDopplerDisplay, SamplesToDisplay
            end
        end

        function [pathDelays, avgPathGains] = interpolate_power_delay_profile_to_system_sample_rate(obj, ...
                                                                                                    pathDelays_beforeInterpolation, ...
                                                                                                    avgPathGains_beforeInterpolation)
            % switch to linear values
            avgPathGains_beforeInterpolation_linear = db2pow(avgPathGains_beforeInterpolation);

            % how often can we sample the power delay profile at our sample rate?
            Ts = 1/obj.config.r_samp_rate;
            n_points = ceil(pathDelays_beforeInterpolation(end)/Ts);
            pathDelays = (0:1:n_points)*Ts;

            % empty container
            avgPathGains_linear = zeros(size(pathDelays));

            % assign each given path delay from the profiles above to the closest sampling point
            for qq = 1:1:numel(pathDelays_beforeInterpolation)
                tmp = pathDelays_beforeInterpolation(qq);
                [~,idx] = min(abs(pathDelays - tmp));
                avgPathGains_linear(idx) = avgPathGains_linear(idx) + avgPathGains_beforeInterpolation_linear(qq);
            end

            % remove any points that have no power as this only increases simulation time to no avail
            tmp_del = (avgPathGains_linear == 0);
            avgPathGains_linear(tmp_del) = []; 
            pathDelays(tmp_del) = [];

            assert(abs(sum(avgPathGains_linear) / sum(avgPathGains_beforeInterpolation_linear) - 1) <= 1e-6);

            % normalize
            avgPathGains_linear = avgPathGains_linear/sqrt(sum(avgPathGains_linear.^2));

            % switch back to dB
            avgPathGains = pow2db(avgPathGains_linear);
        end
        
        function [samples_antenna_rx] = pass_samples(obj, samples_antenna_tx, channel_time_in_seconds)
            
            if nargin == 2
                channel_time_in_seconds = 0;
            end
            
            assert(obj.config.N_TX == size(samples_antenna_tx, 2));
            
            % apply amplitude
            samples_antenna_ch = obj.config.amp*samples_antenna_tx;

            % apply STO, CFO and the carrier phase
            samples_antenna_ch = dectnrp_channel.sto_cfo_phase(samples_antenna_ch, ...
                                                               obj.config.sto_integer, ...
                                                               obj.config.sto_fractional, ...
                                                               obj.config.cfo, ...
                                                               obj.config.err_phase);
                    
            % over-the-air transmission
            if strcmp(obj.config.type, 'AWGN')
                samples_antenna_rx = repmat(sum(samples_antenna_ch, 2), 1, obj.config.N_RX);
            elseif strcmp(obj.config.type, 'Rayleigh') == true || strcmp(obj.config.type, 'Rician') == true
                samples_antenna_rx = obj.r_matlab_MIMO_obj(samples_antenna_ch, channel_time_in_seconds);
            else
                error('Channel neither AWGN, nor Rayleigh, nor Rician, but %s.', obj.config.type);
            end
            
            % thermal noise
            for i=1:1:obj.config.N_RX
                samples_antenna_rx(:,i) = awgn(samples_antenna_rx(:,i), obj.config.snr_db, pow2db(1/obj.config.spectrum_occupied));
            end

            if obj.config.verbosity >= 1
                obj.plot_channel_properties();
            end
        end
        
        function [] = reset_random_Rayleigh_Rician(obj)
            if strcmp(obj.config.type,'Rayleigh') == true || strcmp(obj.config.type,'Rician') == true
                obj.r_matlab_MIMO_obj.reset();
            end
        end

        function [] = plot_channel_properties(obj)
            % lookup PDP
            [pathDelays_beforeInterpolation, avgPathGains_beforeInterpolation] = dectnrp_channel.get_PDP_from_literature(obj.config.r_type, ...
                                                                                                                         obj.config.r_DS_desired);
            avgPathGains_beforeInterpolation_linear = db2pow(avgPathGains_beforeInterpolation);

            % extract for readability
            pathDelays = obj.r_matlab_MIMO_obj.PathDelays;
            avgPathGains = obj.r_matlab_MIMO_obj.AveragePathGains;
            avgPathGains_linear = db2pow(avgPathGains);

            figure()
            clf()

            subplot(2,1,1);
            plot(pathDelays_beforeInterpolation, avgPathGains_beforeInterpolation_linear, 'b-o');
            hold on
            plot(pathDelays, avgPathGains_linear,'r-x');
            title('Channel Path Gains linear');
            xlabel('Time');
            ylabel('Path Gain');
            legend('ITU', 'Interpolation');
            grid on

            % not yet normalized, just for comparison
            avgPathGains = pow2db(avgPathGains_linear);

            subplot(2,1,2);
            plot(pathDelays_beforeInterpolation, avgPathGains_beforeInterpolation, 'b-o');
            hold on
            plot(pathDelays, avgPathGains, 'r-x');
            title('Channel Path Gains logarithmic');
            xlabel('Time');
            ylabel('Path Gain (dB)');
            legend('ITU', 'Interpolation');
            grid on

            % calculate key property for exponential decay
            mean_tau_weighted = sum(pathDelays.*avgPathGains_linear)/sum(avgPathGains_linear);
            rms_tau = sqrt(sum(((pathDelays-mean_tau_weighted).^2).*avgPathGains_linear)/sum(avgPathGains_linear));
        
            offset = 0.74;
            separation = 0.015;
            dectnrp_util.annotation(0.5, offset - 0*separation, sprintf('Sampling Rate: %f MS/s\n', obj.config.r_samp_rate/1e6));
            dectnrp_util.annotation(0.5, offset - 1*separation, sprintf('Sampling time Ts: %f ns\n', 1/obj.config.r_samp_rate/1e-9));
            dectnrp_util.annotation(0.5, offset - 2*separation, sprintf('Largest delay: %f ns\n', pathDelays(end)/1e-9));
            dectnrp_util.annotation(0.5, offset - 3*separation, sprintf('PDP sampling points: %d\n', numel(pathDelays)));
            dectnrp_util.annotation(0.5, offset - 4*separation, sprintf('Tau mean: %f ns\n', mean(pathDelays)/1e-9));
            dectnrp_util.annotation(0.5, offset - 5*separation, sprintf('Tau mean weighted: %f ns\n', mean_tau_weighted/1e-9));
            dectnrp_util.annotation(0.5, offset - 6*separation, sprintf('Tau rms weighted: %f ns\n', rms_tau/1e-9));
            dectnrp_util.annotation(0.5, offset - 7*separation, sprintf('Coherence Bandwidth: %f MHz\n', 1/rms_tau/1e6));
            dectnrp_util.annotation(0.5, offset - 8*separation, sprintf('Occupied Bandwidth: %f MHz\n\n', obj.config.r_samp_rate*obj.config.spectrum_occupied/1e6));
        end
    end
end
