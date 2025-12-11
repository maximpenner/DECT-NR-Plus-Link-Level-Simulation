function [] = run_all()
    clear all;
    close all;
    rng('shuffle');
    warning('on');

    % Tests typically taking < 1 min
    run_single("dectnrp_regression.channel.time_reproducibility");
    run_single("dectnrp_regression.tx_config.channel_coding_pdc");

    % Tests typically taking >> 1 min
    %
    % The number of possible packet configurations is astronomical,
    % which is why it is impossible to test them all. Instead, these
    % tests can either cover a subset or just run endlessly. The 
    % ultimate goal of these regression tests is to trigger an assert,
    % an exception or decode a packet with insufficient SNR, which
    % would indicate a bug in the code.
    run_single("dectnrp_regression.tx_config.single_packet_tx_rx");
    run_single("dectnrp_regression.tx_config.single_packet_tx");

    disp("All tests passed.");
end

function [] = run_single(function_as_str)
    t_Start = tic;
    rng(randi([0 1e9], 1, 1));
    fprintf("Test %s starting with Seed %d.\n", function_as_str, rng().Seed);
    assert(eval(function_as_str));
    fprintf("Test %s finished after %.2f seconds.\n\n", function_as_str, toc(t_Start));
end