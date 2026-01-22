function X = sp_dft(x)
%SP_DFT Direct DFT implementation (O(N^2)) without toolboxes.
%
%   X = sp_dft(x) computes the N-point discrete Fourier transform (DFT)
%   using the definition:
%     X[k] = sum_{n=0}^{N-1} x[n] * exp(-1j*2*pi*k*n/N)
%
%   This matches fft(x) in base MATLAB (up to floating-point roundoff).
%
%   Input:
%     x - Signal (vector, real/complex). Treated as a column vector.
%
%   Output:
%     X - DFT result (N×1 complex column vector)

validateattributes(x, {'numeric'}, {'vector'}, mfilename, 'x', 1);
x = x(:);

N = numel(x);
if N == 0
    X = zeros(0, 1);
    return;
end

x = double(x);
X = zeros(N, 1);

n = (0:N-1).';
for k = 0:(N-1)
    X(k+1) = sum(x .* exp(-1j * 2*pi * (k*n) / N));
end
end

