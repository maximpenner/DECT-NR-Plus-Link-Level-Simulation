%
% Let's say N_b_DFT = 64
%
%   Technical specification (TS) subcarrier indices are:
%
%       31 30 29 .... -26 -27 -28
%
%   Matlab subcarriers are:
%
%       1 2 3 4 5 ... 60
%
%   So this function converts 31 to 1, 30 to 2, -28 to 60 etc.
%   Works in both directions.
%
function [k_i_matlab] = index_conversion_TS_matlab(N_b_DFT, k_i)
    k_i_matlab = N_b_DFT/2 - k_i;
end
