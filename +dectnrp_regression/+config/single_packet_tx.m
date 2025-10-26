function [success] = single_packet_tx()

    range.u_vec                 = [1,2,4,8];
    range.b_vec                 = [1,2,4,8,12,16];
    range.PacketLengthType_vec  = [0,1];
    range.PacketLength_vec      = 2:1:4;
    range.tm_mode_0_to_11_vec   = 0:1:4;
    range.mcs_index_vec         = 1:1:6;
    range.Z_vec                 = [2048,6144];
    range.oversampling_vec      = 1;
    range.codebook_index_vec    = 0:1:0;
    range.PLCF_type_vec         = [1,2];
    range.rv_vec                = 0:1:0;
    range.network_id_vec        = randi([1 2^32], 1, 1);
    range.verbosity             = 0;
    range.inner_repetitions     = 1;

    try
        dectnrp_regression.config.parfor_harness(range, @test_per_config);
    catch ME
        fprintf("Test 0, error message: %s\n\n", dbstack().name, ME.message);
        success = false;
        return;
    end

    success = true;
end

function [] = test_per_config(config)
    tx = dectnrp_tx.tx_t(config);
    tx.generate_random_packet();
end