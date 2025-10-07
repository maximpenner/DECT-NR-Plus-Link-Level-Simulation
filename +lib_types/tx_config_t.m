classdef tx_config_t
    
    properties
        u                   % mu = 1, 2, 4 or 8
        b                   % beta = 1, 2, 4, 8, 12 or 16
        PacketLengthType    % 0 for subslots, 1 for slots
        PacketLength        % min is 1, max is 16 according to Table 6.2.1-2a in part 4
        tm_mode_0_to_11     % Table 7.2-1, mode determines whether transmission is closed loop or not, values range from 0 to 11
        mcs_index           % Table A-1 in part 3, values range from 0 to 11
        Z                   % 5.3, maximum codeblock size
        oversampling        % By how much do we oversample our OFDM packet compared to critical sampling (insert zeros at spectrum edges before IFFT)?
        codebook_index      % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)
        PLCF_type           % Type 1 is 40 bits, Type 2 is 80 bits
        rv                  % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
        network_id          % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY
        verbosity           % show data during execution: 0 false, 1 only text, 2 text + plots
    end
    
    methods
        function obj = tx_config_t()
            obj = obj.set_example_values();
            assert(obj.is_valid());
        end

        function ret = is_valid(obj)
            ret = true;
        end
    end

    methods (Hidden = true)
        function obj = set_example_values(obj)
            obj.u = 1;
            obj.b = 1;
            obj.PacketLengthType = 0;
            obj.PacketLength = 4;
            obj.tm_mode_0_to_11 = 0;
            obj.mcs_index = 5;
            obj.Z = 6144;
            obj.oversampling = 8;
            obj.codebook_index = 0;
            obj.PLCF_type = 2;
            obj.rv = 0;
            obj.network_id = de2bi(1e6,32,'left-msb');
            obj.verbosity = 2;
        end
    end
end
