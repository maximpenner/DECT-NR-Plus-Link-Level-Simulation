function [var_rand] = pm_rand(var)

    sz = size(var);

    % get an evenly distributed random number within -var to +var
    var_rand = ( 2 * (rand(sz)-0.5) ) * var;
end
