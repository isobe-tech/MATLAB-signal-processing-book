function [W, Winv] = sp_dft_matrix(N, varargin)
%SP_DFT_MATRIX Build the DFT matrix W such that X = W*x.
%
%   [W,Winv] = sp_dft_matrix(N) returns:
%     W    - N×N DFT matrix (matches fft scaling)
%     Winv - N×N inverse DFT matrix, so that x = Winv*X
%
%   The default convention matches base MATLAB:
%     X = W*x  corresponds to fft(x)
%     x = (1/N) * W' * X corresponds to ifft(X)
%
%   Name-value options:
%     "Unitary" (default false) - If true, returns a unitary pair:
%       W    = (1/sqrt(N)) * exp(-1j*2*pi*k*n/N)
%       Winv = W'   (so that W'*W = I)

validateattributes(N, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 1);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Unitary", false, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
parse(p, varargin{:});
useUnitary = logical(p.Results.Unitary);

k = (0:N-1).';
n = 0:N-1;
W = exp(-1j * 2*pi * (k*n) / N);

if useUnitary
    W = W / sqrt(N);
    Winv = W';
else
    Winv = (1/N) * W';
end
end

