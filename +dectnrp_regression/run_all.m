function [] = run_all()
    clear all;
    close all;
    clc;

    % test typically taking < 1 min
    run_single("dectnrp_regression.channel.time_reproducibility");
    run_single("dectnrp_regression.config.channel_coding_pdc");

    % test typically taking >> 1 min
    run_single("dectnrp_regression.config.single_packet_tx");
    run_single("dectnrp_regression.config.single_packet_tx_rx");

    disp("All tests passed.");
end

function [] = run_single(function_as_str)
    t_Start = tic;
    rng(randi([0 1e9], 1, 1));
    fprintf("Test %s starting with Seed %d.\n", function_as_str, rng().Seed);
    assert(eval(function_as_str));
    fprintf("Test %s finished after %.2f seconds.\n\n", function_as_str, toc(t_Start));
end