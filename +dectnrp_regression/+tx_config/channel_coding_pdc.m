function [success] = channel_coding_pdc()

    range.u_vec                 = [1,2,4,8];
    range.b_vec                 = [1,2,4,8,12,16];
    range.PacketLengthType_vec  = [0,1];
    range.PacketLength_vec      = 2:1:4;
    range.tm_mode_0_to_11_vec   = 0:1:0;
    range.mcs_index_vec         = 1:1:4;
    range.Z_vec                 = [2048,6144];
    range.oversampling_vec      = 1;
    range.codebook_index_vec    = 0:1:0;
    range.PLCF_type_vec         = 1;
    range.rv_vec                = 0:1:0;
    range.network_id_vec        = randi([1 2^32], 1, 1);
    range.verbosity             = 0;
    range.inner_repetitions     = 1;

    try
        dectnrp_regression.tx_config.parfor_harness(range, @test_per_config);
    catch ME
        s = dbstack;
        fprintf("Test %s, error message: %s\n\n", s(1).name, ME.message);
        success = false;
        return;
    end

    success = true;
end

function [] = test_per_config(tx_config)
    % create transmitter
    tx = dectnrp_tx.tx_t(tx_config);

    % generate random bits
    tb_bits = randi([0 1], tx.tx_derived.N_TB_bits, 1);

    % encode PDC
    [x_PDC, ~] = dectnrp_7_transmission_encoding.PDC_encoding(tb_bits, ...
                                                              tx.tx_derived.G, ...
                                                              tx_config.Z, ...
                                                              tx_config.network_id, ...
                                                              tx_config.PLCF_type, ...
                                                              tx_config.rv, ...
                                                              tx.tx_derived.mcs, ...
                                                              tx.tx_derived.tm_mode.N_SS);

    % simulate without HARQ
    HARQ_buf_ = [];

    % decode PDC
    [tb_bits_recovered, ~, ~] = dectnrp_7_transmission_encoding.PDC_decoding(x_PDC, ...
                                                                             tx.tx_derived.N_TB_bits, ...
                                                                             tx_config.Z, ...
                                                                             tx_config.network_id, ...
                                                                             tx_config.PLCF_type, ...
                                                                             tx_config.rv, ...
                                                                             tx.tx_derived.mcs, ...
                                                                             tx.tx_derived.tm_mode.N_SS, ...
                                                                             HARQ_buf_);

    % analyze result
    n_bit_errors = sum(abs(tb_bits - double(tb_bits_recovered)));

    assert(n_bit_errors == 0);
end