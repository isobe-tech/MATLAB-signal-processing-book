function y = sp_iir1(x, a, varargin)
%SP_IIR1 First-order IIR filter (one-pole) implemented by a difference equation.
%
%   y = sp_iir1(x,a) applies:
%     y[n] = (1-a) x[n] + a y[n-1]
%   with the initial condition y[-1] = 0.
%
%   For 0 < a < 1, this behaves like a simple low-pass smoother:
%   larger a -> slower response (more smoothing, more delay).
%
%   Name-value options:
%     "Y0" (default 0) - Initial output y[-1] (scalar, real/complex)
%
%   Inputs:
%     x - Input signal (vector, real/complex)
%     a - Feedback coefficient (real scalar). Stability requires |a| < 1.
%
%   Output:
%     y - Output signal (N×1 column vector)

validateattributes(x, {'numeric'}, {'vector'}, mfilename, 'x', 1);
validateattributes(a, {'numeric'}, {'scalar', 'real', 'finite'}, mfilename, 'a', 2);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Y0", 0, @(v) isscalar(v) && isnumeric(v) && isfinite(real(v)) && isfinite(imag(v)));
parse(p, varargin{:});
y0 = p.Results.Y0;

x = double(x(:));
N = numel(x);

y = zeros(N, 1);
yp = double(y0);
for n = 1:N
    yc = (1 - a) * x(n) + a * yp;
    y(n) = yc;
    yp = yc;
end
end

