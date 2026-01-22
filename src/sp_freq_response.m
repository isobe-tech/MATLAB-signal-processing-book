function [f, H] = sp_freq_response(fs, b, varargin)
%SP_FREQ_RESPONSE Numerical frequency response evaluation (no toolboxes).
%
%   [f,H] = sp_freq_response(fs,b) evaluates the frequency response of a
%   discrete-time filter with FIR coefficients b at frequencies f [Hz].
%
%   [f,H] = sp_freq_response(fs,b,"a",a,...) evaluates the IIR response
%     H(e^{j\Omega}) = B(e^{j\Omega}) / A(e^{j\Omega})
%   where B and A are polynomials in z^{-1} with coefficients b and a.
%
%   Name-value options:
%     "a"      (default 1)      - Denominator coefficients (scalar or vector)
%     "F"      (default [])     - Frequency points [Hz] (vector)
%     "N"      (default 2000)   - Number of points if "F" is not given
%     "Range"  (default "half") - "half": 0..fs/2, "whole": -fs/2..fs/2
%
%   Outputs:
%     f - Frequency axis [Hz] (column vector)
%     H - Frequency response at f (complex column vector)
%
%   Notes:
%   - This function is intentionally implemented without Signal Processing
%     Toolbox functions such as freqz.
%   - The implementation directly evaluates the DTFT / z-transform on the
%     unit circle using exponentials.
%
%   See also: sp_group_delay

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 1);
b = b(:);
if isempty(b)
    error([mfilename ':EmptyB'], 'b must be non-empty.');
end

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "a", 1, @(v) isnumeric(v) && isvector(v) && ~isempty(v));
addParameter(p, "F", [], @(v) isnumeric(v) && isvector(v));
addParameter(p, "N", 2000, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}));
addParameter(p, "Range", "half", @(v) isstring(v) || ischar(v));
parse(p, varargin{:});

a = p.Results.a(:);
fIn = p.Results.F;
N = p.Results.N;
range = lower(string(p.Results.Range));

if ~isscalar(a) && a(1) == 0
    error([mfilename ':InvalidA'], 'a(1) must be non-zero.');
end

if isempty(fIn)
    switch range
        case "half"
            f = linspace(0, fs/2, N).';
        case "whole"
            f = linspace(-fs/2, fs/2, N).';
        otherwise
            error([mfilename ':InvalidRange'], 'Range must be "half" or "whole".');
    end
else
    f = fIn(:);
end

Omega = 2*pi*f/fs; % rad/sample

mb = 0:numel(b)-1;
B = exp(-1j * (Omega * mb)) * b;

if isscalar(a)
    H = B / a;
else
    ma = 0:numel(a)-1;
    A = exp(-1j * (Omega * ma)) * a;
    H = B ./ A;
end
end

