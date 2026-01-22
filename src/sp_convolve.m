function y = sp_convolve(x, h)
%SP_CONVOLVE Direct convolution (full) with explicit boundary handling.
%
%   y = sp_convolve(x,h) computes the full linear convolution of x and h
%   using the definition (zero-padding outside the given vectors):
%     y[n] = sum_{m=0}^{M-1} h[m] x[n-m]
%   where x has length N, h has length M, and y has length N+M-1.
%
%   This function is written for learning and debugging:
%   - The indexing is explicit (safe for beginners).
%   - The result should match conv(x,h) in base MATLAB.
%
%   Inputs:
%     x - Input signal (vector, real/complex)
%     h - FIR coefficients / impulse response (vector, real/complex)
%
%   Output:
%     y - Full convolution result ((N+M-1)×1 column vector)

validateattributes(x, {'numeric'}, {'vector'}, mfilename, 'x', 1);
validateattributes(h, {'numeric'}, {'vector'}, mfilename, 'h', 2);

x = x(:);
h = h(:);

Nx = numel(x);
Nh = numel(h);

if Nx == 0 || Nh == 0
    y = zeros(0, 1);
    return;
end

% Use double arithmetic for predictable behavior.
x = double(x);
h = double(h);

Ny = Nx + Nh - 1;
y = zeros(Ny, 1);

% y[n] = sum_{m=0..Nh-1} h[m] * x[n-m], with x[k]=0 outside k=0..Nx-1.
for n = 0:(Ny-1)
    acc = 0;
    for m = 0:(Nh-1)
        k = n - m;
        if k >= 0 && k < Nx
            acc = acc + h(m+1) * x(k+1);
        end
    end
    y(n+1) = acc;
end
end

