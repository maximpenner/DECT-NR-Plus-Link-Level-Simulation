function [] = test_tx_packet(tx_config)
    % create transmitter, receiver and channel
    tx = lib_types.tx_t(tx_config);
    
    % create a DECT NR+ packet
    tx.generate_random_packet();
end