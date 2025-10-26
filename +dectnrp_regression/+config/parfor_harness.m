function [] = parfor_harness(range, func)

    assert(isscalar(range.inner_repetitions));
    assert(range.inner_repetitions >= 1);

    % test a large selection of packet configurations
    for u = range.u_vec
        %for b_vec_idx = 1:numel(b_vec)
        parfor b_vec_idx = 1:numel(range.b_vec)
            b = range.b_vec(b_vec_idx);
            for PacketLengthType = range.PacketLengthType_vec
                for PacketLength = range.PacketLength_vec
                    for tm_mode_0_to_11 = range.tm_mode_0_to_11_vec
                        for mcs_index = range.mcs_index_vec
                            for Z = range.Z_vec
                                for oversampling = range.oversampling_vec
                                    for codebook_index = range.codebook_index_vec
                                        for PLCF_type = range.PLCF_type_vec
                                            for rv = range.rv_vec
                                                for network_id = range.network_id_vec
                                                    for inner_repetition = 1:range.inner_repetitions
    
                                                        % configurations are initialized with exemplary values
                                                        config = dectnrp_tx.config_t();
    
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
                                                        config.verbosity = range.verbosity;
                            
                                                        func(config);
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
end