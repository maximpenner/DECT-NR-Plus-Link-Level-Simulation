function [b, segErr] = Code_block_desegmentation_Z_2048(c_r, B)

    % source: https://github.com/robmaunder/turbo-3gpp-matlab
    
    % according to 5.3 filler bits are unnecessary

    assert(B > 0, 'Unsupported block length');

    supported_values_of_K = [40:8:511,512:16:1023,1024:32:2047,2048:64:6144];

    Z = 2048;

    if B <= Z
        L = 0;
        C = 1;
        B_prime = B;
    else
        L = 24;
        C = ceil(B/(Z-L));
        B_prime = B+C*L;
    end

    K_plus = min(supported_values_of_K(C*supported_values_of_K>=B_prime));

    if C == 1
        C_plus = 1;
        K_minus = 0;
        C_minus = 0;
    elseif C>1
        K_minus = max(supported_values_of_K(supported_values_of_K<K_plus));
        delta_K = K_plus - K_minus;
        C_minus = floor((C*K_plus-B_prime)/delta_K);
        C_plus = C-C_minus;
    end

    K_r = zeros(1,C);
    for r = 0:C-1
        if r < C_minus
            K_r(r+1) = K_minus;
        else
            K_r(r+1) = K_plus;
        end
    end

    F = C_plus*K_plus + C_minus*K_minus - B_prime;
    
    assert(F == 0);

    b = zeros(1,B);
    segErr = [];

    k = F;
    s = 0;
    for r = 0:C-1
        while k < K_r(r+1)-L
            b(s+1) = c_r{r+1}(k+1);
            k = k+1;
            s = s+1;
        end

        if C>1
            a_r = c_r{r+1}(1:K_r(r+1)-L);
            
            assert(sum(isnan(a_r)) == 0);
            
            a_r(isnan(a_r)) = 0;

            %p_r2 = calculate_crc_bits(a_r,G_max);
            p_r2 = lteCRCEncode(a_r,'24B');
            p_r2 = p_r2(end-23:end);
            p_r2 = double(p_r2');

            p_r = zeros(1,L);
            while k < K_r(r+1)
                p_r(k+L-K_r(r+1)+1) = c_r{r+1}(k+1);
                k = k+1;
            end

            if ~isequal(p_r,p_r2)
                %b = [];
                segErr = [segErr int8(1)];
                %return;
            end

        end
        k=0;
    end
    
    % convert b to correct format
    b = b';
    b = int8(b);    
end
