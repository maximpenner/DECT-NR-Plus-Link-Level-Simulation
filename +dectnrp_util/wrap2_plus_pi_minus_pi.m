function [angles] = wrap2_plus_pi_minus_pi(angles)
    % https://stackoverflow.com/questions/28313558/how-to-wrap-a-number-into-a-range
    angles = mod(angles + pi, 2*pi) - pi;
end
