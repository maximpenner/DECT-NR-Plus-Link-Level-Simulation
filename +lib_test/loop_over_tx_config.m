function [] = loop_over_tx_config()

    % Note that some errors are to be expected. For example,
    % there are different numbers of codebook indices for the
    % antenna configurations, and some packet lengths for 1
    % or 2 antennas are too short when using a configuration
    % with 4 or 8 antennas.

    % configurations are initialized with exemplary values
    tx_config = lib_types.tx_config_t();

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

    cnt = 0;

    % test a large selection of packet configurations
    for u = u_vec
        for b = b_vec
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

                                                    % overwrite
                                                    tx_config.u = u;
                                                    tx_config.b = b;
                                                    tx_config.PacketLengthType = PacketLengthType;
                                                    tx_config.PacketLength = PacketLength;
                                                    tx_config.tm_mode_0_to_11 = tm_mode_0_to_11;
                                                    tx_config.mcs_index = mcs_index;
                                                    tx_config.Z = Z;
                                                    tx_config.oversampling = oversampling;
                                                    tx_config.codebook_index = codebook_index;
                                                    tx_config.PLCF_type = PLCF_type;
                                                    tx_config.rv = rv;
                                                    tx_config.network_id = de2bi(network_id,32,'left-msb');
                                                    tx_config.verbosity = verbosity;
                                                    
                                                    % test channel coding
                                                    try
                                                        result = lib_test.test_channel_coding_pdc(tx_config);
                                                        assert(result.n_bit_errors == 0)
                                                    catch ME
                                                        fprintf("Error for cnt=%d. Message: %s\n", cnt, ME.message);
                                                        return;
                                                    end
                        
                                                    % test packet generation
                                                    try
                                                        lib_test.test_tx_packet(tx_config);
                                                    catch ME
                                                        fprintf("Error for cnt=%d. Message: %s\n", cnt, ME.message);
                                                        return;
                                                    end

                                                    cnt = cnt + 1;
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
