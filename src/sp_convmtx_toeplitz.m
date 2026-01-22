function H = sp_convmtx_toeplitz(h, Nx)
%SP_CONVMTX_TOEPLITZ Convolution matrix (Toeplitz) for full FIR convolution.
%
%   H = sp_convmtx_toeplitz(h,Nx) returns a matrix H of size
%   (Nx+M-1)×Nx (M = numel(h)) such that:
%     y = H * x
%   matches the full convolution y = conv(x,h) for any x of length Nx.
%
%   This is a toolbox-free alternative to convmtx (DSP/Signal Processing TB).
%   It is intended for understanding the "filter = linear transform" view.
%
%   Inputs:
%     h  - FIR coefficients / impulse response (vector, length M)
%     Nx - Input length (positive integer scalar)
%
%   Output:
%     H - Convolution matrix ((Nx+M-1)×Nx)

validateattributes(h, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'h', 1);
validateattributes(Nx, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'Nx', 2);

h = double(h(:));
M = numel(h);

H = zeros(Nx + M - 1, Nx);
for i = 1:Nx
    H(i:(i+M-1), i) = h;
end
end

