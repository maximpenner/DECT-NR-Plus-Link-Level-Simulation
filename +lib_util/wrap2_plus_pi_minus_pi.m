function [angles] = wrap2_plus_pi_minus_pi(angles)

    % https://stackoverflow.com/questions/28313558/how-to-wrap-a-number-into-a-range

    % Since % operator lets you to get a [0; limit) space, a standard approach when you need a [start; limit) space is to get a space [0; limit - start] via % operator and then add it to start.
    % This way you'll get a [start; limit) space. So the formula would be
    % 
    % numberToBeWrapped = start + (numberToBeWrapped - start) % (limit - start)

    angles = mod(angles + pi, 2*pi) - pi;
end
