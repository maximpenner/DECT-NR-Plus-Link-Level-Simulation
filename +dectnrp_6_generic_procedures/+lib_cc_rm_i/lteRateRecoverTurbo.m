function d_turbo = lteRateRecoverTurbo(f, N_TB_bits, rv, chs, cbsbuffers, G, N_L, Q_m)
    % The function lteRateRecoverTurbo() provided by Matlab is given the transport block size as an argument.
    % Based on the transport block size, it deduces the code block sizes assuming Z=6144.
    % To allow Z=2048, this same-named function processes each code block individually.
    % This way code block sizes for Z=2048 can be processed.

    % sizes of all code blocks before and after encoding
    K_r = dectnrp_6_generic_procedures.lib_cc_rm_i.get_3gpp_code_block_segment_lengths(N_TB_bits + 24, 2048);
    C = numel(K_r);
    E_r = dectnrp_6_generic_procedures.lib_cc_rm_i.get_3gpp_encoded_code_block_segment_lengths(G, C, N_L, Q_m);

    assert(numel(K_r) > 0);
    assert(numel(K_r) == C);
    assert(numel(K_r) == numel(E_r));
    assert(numel(f) == sum(E_r));

    d_turbo = cell(1, numel(E_r));

    % process each code block individually
    f_cnt = 1;
    for i=1:1:numel(E_r)

        % extract sub range of bits
        f_sub = f(f_cnt : f_cnt + E_r(i) - 1);
        f_cnt = f_cnt + E_r(i);

        assert(numel(f_sub) == E_r(i));

        % This is the transport block size that is given lteRateRecoverTurbo() as an argument.
        % The function lteRateRecoverTurbo() will assume there is only one code block and add the 24 bit CRC.
        trblklen = K_r(i) - 24;

        if numel(cbsbuffers) > 0
            d_turbo(i) = lteRateRecoverTurbo(f_sub, trblklen, rv, chs, cbsbuffers(i));
        else
            d_turbo(i) = lteRateRecoverTurbo(f_sub, trblklen, rv, chs);
        end
    end
end