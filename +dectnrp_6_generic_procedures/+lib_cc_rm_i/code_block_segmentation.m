function [c_r] = code_block_segmentation(b)

    % source: https://github.com/robmaunder/turbo-3gpp-matlab
    
    % according to 5.3 filler bits are unnecessary
    
    B = numel(b);

    assert(B > 0, 'Unsupported block length');

    supported_values_of_K = [40:8:511,512:16:1023,1024:32:2047,2048];

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

    c_r = cell(1,C);
    for r = 0:C-1
        c_r{r+1} = zeros(1,K_r(r+1));
    end

    for k = 0:F-1
        c_r{1}(k+1) = NaN;
    end

    k = F;
    s = 0;
    for r = 0:C-1
        while k < K_r(r+1)-L
            c_r{r+1}(k+1) = b(s+1);
            k = k+1;
            s = s+1;
        end

        if C>1
            a_r = c_r{r+1}(1:K_r(r+1)-L);
            
            assert(sum(isnan(a_r)) == 0);
            
            a_r(isnan(a_r)) = 0;

            %p_r = calculate_crc_bits(a_r,G_max);
            temp = lteCRCEncode(a_r,'24B');
            p_r = temp(end-23:end);
            p_r = p_r';

            while k < K_r(r+1)
                c_r{r+1}(k+1) = p_r(k+L-K_r(r+1)+1);
                k = k+1;
            end
        end
        k=0;
        
        % we need to convert to the same format as matlab uses
        temp = cell2mat(c_r(r+1));
        temp = temp';
        temp = int8(temp);
        c_r(r+1) = {temp};
    end
end
