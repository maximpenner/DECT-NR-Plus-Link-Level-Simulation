function [] = plot_STF(samples_antenna_tx, u, N_b_DFT, oversampling)

    N_b_DFT = N_b_DFT*oversampling;

    % extract length of repetitive pattern
    switch u
        case 1
            STF_len = N_b_DFT*9/8*14/9;
            rep_seq = 7;
        case {2,4,8}
            STF_len = N_b_DFT*9/8*2;
            rep_seq = 9;
    end
    seq_len = STF_len/rep_seq;
    
    % abs
    figure()
    clf()

    subplot(3,2,1)
    for i=1:rep_seq
        seq = abs( samples_antenna_tx((i-1)*seq_len+1:i*seq_len,1) );
        plot(abs(seq));
        hold on;
    end
    ylabel('Pattern abs')
    legend();
    
    % real
    subplot(3,2,3)
    for i=1:rep_seq
        seq = real( samples_antenna_tx((i-1)*seq_len+1:i*seq_len,1) );
        plot(abs(seq));
        hold on;
    end
    ylabel('Pattern real')
    
    % imag
    subplot(3,2,5)
    for i=1:rep_seq
        seq = imag( samples_antenna_tx((i-1)*seq_len+1:i*seq_len,1) );
        plot(abs(seq));
        hold on;
    end
    ylabel('Pattern imag')
    
    % plot a few more sample to show that the pattern does not repeat further
    i = rep_seq + 1;
    subplot(3,2,2)
    seq = abs( samples_antenna_tx((i-1)*seq_len+1:i*seq_len,1) );
    plot(abs(seq));
    ylabel('next symbols abs')
    
    subplot(3,2,4)
    seq = real( samples_antenna_tx((i-1)*seq_len+1:i*seq_len,1) );
    plot(abs(seq));
    ylabel('next symbols real')
    
    subplot(3,2,6)
    seq = imag( samples_antenna_tx((i-1)*seq_len+1:i*seq_len,1) );
    plot(abs(seq));
    ylabel('next symbols imag')

    sgtitle('STF repetitive pattern');
end

