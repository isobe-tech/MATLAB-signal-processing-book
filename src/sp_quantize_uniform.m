function [xq, qerr, step, code] = sp_quantize_uniform(x, bits, varargin)
%SP_QUANTIZE_UNIFORM Uniform quantizer (round-to-nearest).
%
%   xq = sp_quantize_uniform(x,bits) quantizes x to a uniform grid whose
%   step is determined by bits and FullScale (default 1). Quantization is
%   performed by rounding to the nearest step and saturating to the
%   representable range.
%
%   [xq,qerr,step,code] = ... also returns quantization error qerr = xq - x,
%   the step size, and integer code indices.
%
%   Name-value options:
%     "FullScale" (default 1)     - Full-scale (peak) value (> 0)
%     "Saturate"  (default true)  - Saturate to [-FullScale, +FullScale]
%
%   Notes:
%     This model intentionally prioritizes simplicity for learning. The exact
%     endpoint behavior of practical ADC/DAC quantizers may differ by 1 LSB.

validateattributes(bits, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'bits', 2);
x = x(:);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "FullScale", 1, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
addParameter(p, "Saturate", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
parse(p, varargin{:});

fullScale = double(p.Results.FullScale);
saturate = logical(p.Results.Saturate);

% Step corresponding to B bits over a peak-to-peak range of 2*FullScale.
step = (2 * fullScale) / (2^double(bits));

code = round(x / step);
codeMin = -2^(double(bits) - 1);
codeMax = 2^(double(bits) - 1);

if saturate
    code = min(max(code, codeMin), codeMax);
end

xq = code * step;
qerr = xq - x;
end

