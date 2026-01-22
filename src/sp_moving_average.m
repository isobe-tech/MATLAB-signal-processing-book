function [y, h] = sp_moving_average(x, M)
%SP_MOVING_AVERAGE Causal moving average FIR filter (length-N output).
%
%   y = sp_moving_average(x,M) applies a causal moving average with window
%   length M (zero initial condition). The output y has the same length as x.
%
%     h[m] = 1/M,  m = 0..M-1
%     y[n] = sum_{m=0..M-1} h[m] x[n-m], with x[k]=0 for k<0
%
%   [y,h] = ... also returns the FIR coefficients h (M×1).
%
%   Input:
%     x - Input signal (vector)
%     M - Window length (positive integer scalar)
%
%   Output:
%     y - Output signal (N×1)
%     h - FIR coefficients (M×1)

validateattributes(x, {'numeric'}, {'vector'}, mfilename, 'x', 1);
validateattributes(M, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'M', 2);

x = double(x(:));
N = numel(x);

h = ones(M, 1) / M;
yFull = sp_convolve(x, h);
y = yFull(1:N);
end

