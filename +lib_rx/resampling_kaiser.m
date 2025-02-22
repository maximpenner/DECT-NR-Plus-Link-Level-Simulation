function fir_coef = resampling_kaiser(  PassbandFrequency,...
                                        StopbandFrequency,...
                                        PassbandRipple_dB,...
                                        StopbandAttenuation_dB,...
                                        SampleRate,...
                                        N_force_odd)

    % sources
    %   https://de.mathworks.com/help/signal/ref/kaiser.html
    %   https://de.mathworks.com/help/signal/ref/kaiserord.html
    %   https://tomroelandts.com/articles/how-to-create-a-simple-low-pass-filter
    %   https://tomroelandts.com/articles/how-to-create-a-configurable-filter-using-a-kaiser-window

    % normalize frequencies
    PassbandFrequency = PassbandFrequency/SampleRate;
    StopbandFrequency = StopbandFrequency/SampleRate;
    
    % different terminology in sources
    b           = StopbandFrequency - PassbandFrequency;    % transition bandwidth
    f_cutoff    = PassbandFrequency + b/2;                  % cutoff frequency in the center between passband and stopband
    
    % The filter coefficients are designed by using A, which is the stopband attenuation.
    % The factor A can defined by using StopbandAttenuation_dB, or it can be derived from PassbandRipple_dB depending of which requirement is stricter (the smaller the passband ripple, the larger the stopband attenuation).
    % Thereore, we have to check which is is the strictest definition of A (the largest value of A).
    delta_option_0 = 10^(-StopbandAttenuation_dB/20);
    delta_option_1 = 10^(PassbandRipple_dB/20)-1;
    delta = min(delta_option_0, delta_option_1);
    A = -20*log10(delta);
    
    % Kaiser window parameter beta derived from A
    if A > 50
        beta = 0.1102*(A-8.7);
    elseif 21 <= A && A <= 50
        beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
    elseif A < 21
        beta = 0;
    end
        
    % what filter order is required to fullfill requirements?
    filter_order = (A-7.95)/(2.285*2*pi*b);

    % number of taps
    N = ceil(filter_order + 1);

    assert(N > 0);

    % make sure we have an uneven filter length, this way the delay isn't fractional
    if N_force_odd == true
        if mod(N,2) == 0
            N = N + 1;
        end
    end
    
    % time base
    n = 0:1:(N-1);
    
    % Kaiser window coefficients
    w_n = beta*sqrt(1-(2*n/(N-1)-1).^2);

    % matlab function and custom function
    %w_n = besseli(0,w_n)/besseli(0,beta);
    tmp0 = modified_bessel_function_I_0(beta);
    for i=1:1:N
        w_n(i) = modified_bessel_function_I_0(w_n(i))/tmp0;
    end
    
    % sinc filter centered at 0
    h_n = 2*f_cutoff*sinc(2*f_cutoff*(n-(N-1)/2));
    
    % final coefficients
    fir_coef = w_n.*h_n;

    % normalize
    fir_coef = fir_coef/sum(fir_coef);
end

% source: last equation in https://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html
function val = modified_bessel_function_I_0(z)

    % for implementation in C++ use double, convert final coefficients to float

    % 10 is enough for very precise coefficients
    k_max = 8;

    val = zeros(k_max,1);

    for k=0:1:k_max

        A = (0.25*z^2)^k;

        % eqauivalent to k!
        B = 1;
        for i=2:1:k
            B = B*i;
        end

        B = B*B;

        val(k+1) = A/B;
    end

    val = sum(val);
end