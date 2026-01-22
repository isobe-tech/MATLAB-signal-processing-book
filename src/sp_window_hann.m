function w = sp_window_hann(N)
%SP_WINDOW_HANN Hann window, as a column vector (no toolboxes).
%
%   w = sp_window_hann(N) returns an N-by-1 Hann window:
%     w[n] = 0.5 - 0.5*cos(2*pi*n/(N-1)),  n = 0..N-1
%
%   Notes:
%   - For N=1, the window is defined as w = 1.

validateattributes(N, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 1);

if N == 1
    w = 1;
    return;
end

n = (0:N-1).';
w = 0.5 - 0.5*cos(2*pi*n/(N-1));
end





