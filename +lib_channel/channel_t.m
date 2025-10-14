classdef channel_t < handle
    
    properties
        
        verbosity;
        type;                   % awgn, rayleigh, rician
        
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

        % prefix r_ for rayleigh and rician
        r_samp_rate;            % system sample rate in Samples/s
        r_max_doppler;          % maximum doppler in Hertz
        r_type;                 % Exponential Decay: TDL-i, TDL-ii, TDL-iii(NLOS) or TDL-iv, TDL-v (LOS) etc.
        r_DS_desired;           % scaling factor of normalized delay spread (e.g. according to ITU-R M.2412-0, Table A1-43)
        r_K;                    % rician fading K factor
        
        % variables set and used internally
        r_matlab_MIMO_obj;      % reference to the matlab object -> will be set in init function
        r_appendix;             % appended samples by rayleigh channel -> will be set in init function
    end
    
    methods (Static = true, Access = public)
       
       % Init parameters outside of constructor, or use the rf_channel_example_factory.
       % This is typically more convenient than controlling all parameters through class methods.
       function obj = channel_t()
            obj.verbosity               = 0;
            obj.type                    = [];
            
            obj.N_TX                    = [];
            obj.N_RX                    = [];
            obj.spectrum_occupied       = [];

            obj.amp                     = [];
            obj.sto_integer             = [];
            obj.sto_fractional          = [];
            obj.cfo                     = [];
            obj.err_phase               = [];
            obj.snr_db                  = [];
            
            obj.r_samp_rate             = [];
            obj.r_max_doppler           = [];
            obj.r_type                  = [];
            obj.r_DS_desired            = [];
            obj.r_K                     = [];          

            obj.r_matlab_MIMO_obj       = [];
            obj.r_appendix              = [];
       end
    end
    
    methods (Static = false, Access = public)
        function init_rayleigh_rician_channel(obj)
            
            if strcmp(obj.type,'rayleigh') == false && strcmp(obj.type,'rician') == false
                if strcmp(obj.type,'awgn')
                    return;
                else
                    error('Unknown channel type %s. Must be awgn, rayleigh or rician.', obj.type);
                end              
            end
            
            % lookup PDP
            [pathDelays_beforeInterpolation, avgPathGains_beforeInterpolation] = lib_channel.get_PDP_from_literature(obj.r_type, ...
                                                                                                                     obj.r_DS_desired);
            
            % interpolate to system sample rate
            [pathDelays, avgPathGains] = obj.interpolate_power_delay_profile_to_system_sample_rate(pathDelays_beforeInterpolation, ...
                                                                                                   avgPathGains_beforeInterpolation);
            
            assert(numel(pathDelays) == numel(unique(pathDelays)));
            assert(sum(pathDelays < 0) == 0);
            assert(sum(isnan(avgPathGains)) == 0);
            assert(sum(~isfinite(avgPathGains)) == 0);
            
            % determine how many samples the channel is adding
            Ts = 1/obj.r_samp_rate;
            obj.r_appendix = ceil(max(pathDelays)/Ts);

            % create Matlab channel object
            if strcmp(obj.type,'rayleigh') == true
                
                assert(strcmp(obj.r_type, 'TDL-i') || ...
                       strcmp(obj.r_type, 'TDL-ii') || ...
                       strcmp(obj.r_type, 'TDL-iii') || ...
                       strcmp(obj.r_type, 'InH NLOS') || ...
                       strcmp(obj.r_type, 'UMi NLOS') || ...
                       strcmp(obj.r_type, 'UMi OtI'));

                obj.r_matlab_MIMO_obj = comm.MIMOChannel('SampleRate', obj.r_samp_rate, ...
                                                         'PathDelays', pathDelays, ...
                                                         'AveragePathGains', avgPathGains, ...
                                                         'NormalizePathGains', true, ...
                                                         'FadingDistribution', 'Rayleigh', ... % KFactor, DirectPathDopplerShift, DirectPathInitialPhase
                                                         'MaximumDopplerShift', obj.r_max_doppler, ...
                                                         'DopplerSpectrum', doppler('Jakes'), ...
                                                         'SpatialCorrelationSpecification', 'Separate Tx Rx', ... %'NumTransmitAntennas', obj.N_TX, 'NumReceiveAntennas', obj.N_RX, ...
                                                         'TransmitCorrelationMatrix', eye(obj.N_TX), ...
                                                         'ReceiveCorrelationMatrix', eye(obj.N_RX), ... % SpatialCorrelationMatrix
                                                         'AntennaSelection', 'off', ...
                                                         'NormalizeChannelOutputs', false, ...
                                                         'FadingTechnique', 'Sum of sinusoids', ...
                                                         'NumSinusoids', 48, ...
                                                         'InitialTimeSource', 'Input Port', ... % InitialTime
                                                         'RandomStream', 'Global stream', ... % Seed
                                                         'PathGainsOutputPort', false, ...
                                                         'Visualization', 'off'); %AntennaPairsToDisplay, PathsForDopplerDisplay, SamplesToDisplay
                
            elseif strcmp(obj.type,'rician') == true
                
                assert(strcmp(obj.r_type, 'TDL-iv') || ...
                       strcmp(obj.r_type, 'TDL-v') || ...
                       strcmp(obj.r_type, 'InH LOS') || ...
                       strcmp(obj.r_type, 'UMi LOS'));
                
                obj.r_matlab_MIMO_obj = comm.MIMOChannel('SampleRate', obj.r_samp_rate, ...
                                                         'PathDelays', pathDelays, ...
                                                         'AveragePathGains', avgPathGains, ...
                                                         'NormalizePathGains', true, ...
                                                         'FadingDistribution', 'Rician', ...
                                                         'KFactor', obj.r_K, ...
                                                         'DirectPathDopplerShift',obj.r_max_doppler*0.7, ...
                                                         'DirectPathInitialPhase',0, ...
                                                         'MaximumDopplerShift', obj.r_max_doppler, ...
                                                         'DopplerSpectrum', doppler('Jakes'), ...
                                                         'SpatialCorrelationSpecification', 'Separate Tx Rx', ... %'NumTransmitAntennas', obj.N_TX, 'NumReceiveAntennas', obj.N_RX, ...
                                                         'TransmitCorrelationMatrix', eye(obj.N_TX), ...
                                                         'ReceiveCorrelationMatrix', eye(obj.N_RX), ... % SpatialCorrelationMatrix
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
            Ts = 1/obj.r_samp_rate;
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
            
            assert(obj.N_TX == size(samples_antenna_tx, 2));
            
            % apply amplitude
            samples_antenna_ch = obj.amp*samples_antenna_tx;

            % apply STO, CFO and the carrier phase
            samples_antenna_ch = lib_channel.sto_cfo_phase(samples_antenna_ch, ...
                                                           obj.sto_integer, ...
                                                           obj.sto_fractional, ...
                                                           obj.cfo, ...
                                                           obj.err_phase);
            
            % over-the-air transmission
            if strcmp(obj.type, 'awgn')
                samples_antenna_rx = repmat(sum(samples_antenna_ch, 2), 1, obj.N_RX);
            elseif strcmp(obj.type, 'rayleigh') == true || strcmp(obj.type, 'rician') == true
                samples_antenna_rx = obj.r_matlab_MIMO_obj(samples_antenna_ch, channel_time_in_seconds);
            else
                error('Channel neither awgn, nor rayleigh, nor rician, but %s.', obj.type);
            end
            
            % thermal noise
            for i=1:1:obj.N_RX
                samples_antenna_rx(:,i) = awgn(samples_antenna_rx(:,i), obj.snr_db, pow2db(1/obj.spectrum_occupied));
            end

            if obj.verbosity >= 1
                obj.plot_channel_properties();
            end
        end
        
        function [] = reset_random_rayleigh_rician(obj)
            if strcmp(obj.type,'rayleigh') == true || strcmp(obj.type,'rician') == true
                obj.r_matlab_MIMO_obj.reset();
            end
        end

        function [] = plot_channel_properties(obj)
            % lookup PDP
            [pathDelays_beforeInterpolation, avgPathGains_beforeInterpolation] = lib_channel.get_PDP_from_literature(obj.r_type, ...
                                                                                                                     obj.r_DS_desired);
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
            lib_util.annotation(0.5, offset - 0*separation, sprintf('Sampling Rate: %f MS/s\n', obj.r_samp_rate/1e6));
            lib_util.annotation(0.5, offset - 1*separation, sprintf('Sampling time Ts: %f ns\n', 1/obj.r_samp_rate/1e-9));
            lib_util.annotation(0.5, offset - 2*separation, sprintf('Largest delay: %f ns\n', pathDelays(end)/1e-9));
            lib_util.annotation(0.5, offset - 3*separation, sprintf('PDP sampling points: %d\n', numel(pathDelays)));
            lib_util.annotation(0.5, offset - 4*separation, sprintf('Tau mean: %f ns\n', mean(pathDelays)/1e-9));
            lib_util.annotation(0.5, offset - 5*separation, sprintf('Tau mean weighted: %f ns\n', mean_tau_weighted/1e-9));
            lib_util.annotation(0.5, offset - 6*separation, sprintf('Tau rms weighted: %f ns\n', rms_tau/1e-9));
            lib_util.annotation(0.5, offset - 7*separation, sprintf('Coherence Bandwidth: %f MHz\n', 1/rms_tau/1e6));
            lib_util.annotation(0.5, offset - 8*separation, sprintf('Occupied Bandwidth: %f MHz\n\n', obj.r_samp_rate*obj.spectrum_occupied/1e6));
        end
    end
end
