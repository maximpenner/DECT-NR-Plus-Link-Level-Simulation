function [] = loop_over_config()
    % Note that some errors are to be expected. For example,
    % there are different numbers of codebook indices for the
    % antenna configurations, and some packet lengths for 1
    % or 2 antennas are too short when using a configuration
    % with 4 or 8 antennas.

    % test range
    u_vec                   = [1,2,4,8];
    b_vec                   = [1,2,4,8,12,16];
    PacketLengthType_vec    = [0,1];
    PacketLength_vec        = 2:1:4;
    tm_mode_0_to_11_vec     = 0:1:0;
    mcs_index_vec           = 1:1:4;
    Z_vec                   = [2048,6144];
    oversampling_vec        = [1,2,4,8];
    codebook_index_vec      = 0:1:1;
    PLCF_type_vec           = [1,2];
    rv_vec                  = 0:1:0;
    network_id_vec          = randi([1 2^32], 1, 1);
    verbosity               = 0;

    % test a large selection of packet configurations
    for u = u_vec
        %for b_vec_idx = 1:numel(b_vec)
        parfor b_vec_idx = 1:numel(b_vec)
            b = b_vec(b_vec_idx);
            for PacketLengthType = PacketLengthType_vec
                for PacketLength = PacketLength_vec
                    for tm_mode_0_to_11 = tm_mode_0_to_11_vec
                        for mcs_index = mcs_index_vec
                            for Z = Z_vec
                                for oversampling = oversampling_vec
                                    for codebook_index = codebook_index_vec
                                        for PLCF_type = PLCF_type_vec
                                            for rv = rv_vec
                                                for network_id = network_id_vec

                                                    % configurations are initialized with exemplary values
                                                    config = lib_tx.config_t();

                                                    % overwrite
                                                    config.u = u;
                                                    config.b = b;
                                                    config.PacketLengthType = PacketLengthType;
                                                    config.PacketLength = PacketLength;
                                                    config.tm_mode_0_to_11 = tm_mode_0_to_11;
                                                    config.mcs_index = mcs_index;
                                                    config.Z = Z;
                                                    config.oversampling = oversampling;
                                                    config.codebook_index = codebook_index;
                                                    config.PLCF_type = PLCF_type;
                                                    config.rv = rv;
                                                    config.network_id = de2bi(network_id,32,'left-msb');
                                                    config.verbosity = verbosity;
                                                    
                                                    % test channel coding
                                                    try
                                                        result = lib_regression.test_channel_coding_pdc(config);
                                                        assert(result.n_bit_errors == 0)
                                                    catch ME
                                                        disp(config);
                                                        fprintf("Error for cnt=%d. Message: %s\n\n", ME.message);
                                                    end
                        
                                                    % test packet generation
                                                    try
                                                        lib_regression.test_single_packet(config);
                                                    catch ME
                                                        disp(config);
                                                        fprintf("Error for cnt=%d. Message: %s\n\n", ME.message);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    fprintf("No errors for cnt=%d.\n", cnt);
end
