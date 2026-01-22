function w = sp_window_hamming(N)
%SP_WINDOW_HAMMING Hamming window, as a column vector (no toolboxes).
%
%   w = sp_window_hamming(N) returns an N-by-1 Hamming window:
%     w[n] = 0.54 - 0.46*cos(2*pi*n/(N-1)),  n = 0..N-1
%
%   Notes:
%   - For N=1, the window is defined as w = 1.

validateattributes(N, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 1);

if N == 1
    w = 1;
    return;
end

n = (0:N-1).';
w = 0.54 - 0.46*cos(2*pi*n/(N-1));
end





