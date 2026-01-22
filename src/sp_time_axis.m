function [t, n] = sp_time_axis(fs, N)
%SP_TIME_AXIS Time axis and sample indices for discrete-time signals.
%
%   [t,n] = sp_time_axis(fs,N) returns column vectors (N×1):
%     n = (0:N-1)'
%     t = n/fs  [s]
%
%   Inputs:
%     fs - Sampling frequency [Hz] (positive scalar)
%     N  - Number of samples (positive integer scalar)
%
%   Outputs:
%     t - Time axis [s] (N×1 column vector)
%     n - Sample indices (N×1 column vector)

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 1);
validateattributes(N, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 2);

n = (0:N-1).';
t = n / fs;
end
