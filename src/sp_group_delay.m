function [f, gd] = sp_group_delay(fs, b, varargin)
%SP_GROUP_DELAY Approximate group delay from the frequency response.
%
%   [f,gd] = sp_group_delay(fs,b) returns the approximate group delay gd
%   [samples] of an FIR filter b, evaluated on a frequency grid f [Hz].
%
%   [f,gd] = sp_group_delay(fs,b,"a",a,...) supports an IIR denominator a
%   (scalar or vector) and uses the same frequency grid options as
%   sp_freq_response.
%
%   The group delay is defined as:
%     gd(\Omega) = - d/d\Omega arg(H(e^{j\Omega}))
%   This implementation uses phase unwrapping followed by a numerical
%   derivative.
%
%   Name-value options:
%     "a"       (default 1)       - Denominator coefficients (scalar or vector)
%     "F"       (default [])      - Frequency points [Hz] (vector)
%     "N"       (default 2000)    - Number of points if "F" is not given
%     "Range"   (default "half")  - "half": 0..fs/2, "whole": -fs/2..fs/2
%     "MinMag"  (default 1e-6)    - Mask gd where |H| is too small
%
%   Outputs:
%     f  - Frequency axis [Hz] (column vector)
%     gd - Group delay [samples] (column vector; NaN where undefined)
%
%   See also: sp_freq_response

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
addParameter(p, "MinMag", 1e-6, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'}));
parse(p, varargin{:});

a = p.Results.a;
minMag = p.Results.MinMag;

[f, H] = sp_freq_response(fs, b, "a", a, "F", p.Results.F, "N", p.Results.N, "Range", p.Results.Range);

Omega = 2*pi*f/fs;
phi = unwrap(angle(H));
gd = -gradient(phi, Omega); % samples

mag = abs(H);
gd(mag < minMag) = NaN;
end

