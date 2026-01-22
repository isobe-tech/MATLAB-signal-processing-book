function spec = sp_filter_specs_eval(f, H, varargin)
%SP_FILTER_SPECS_EVAL Evaluate basic filter specs from frequency response (no toolboxes).
%
%   spec = sp_filter_specs_eval(f,H,...) evaluates passband ripple and stopband
%   attenuation from a complex frequency response H sampled on frequency axis f [Hz].
%
%   Name-value options:
%     "Passband"      (default []) - Passband interval(s) [Hz], K×2 matrix [f1 f2]
%     "Stopband"      (default []) - Stopband interval(s) [Hz], K×2 matrix [f1 f2]
%     "PassbandGain"  (default 1)  - Reference gain for passband ripple
%     "Eps"           (default 1e-20) - Added before log to avoid -Inf
%
%   Output struct fields (when the corresponding band is provided):
%     passbandMaxDb, passbandMinDb, passbandRippleDb
%     stopbandMaxDb, stopbandAttenDb
%
%   Notes:
%   - Magnitude (amplitude) specs use 20*log10 (not 10*log10).
%   - This function is intentionally simple: "how to measure" is part of learning.
%
%   See also: sp_freq_response

f = f(:);
H = H(:);
if numel(f) ~= numel(H)
    error([mfilename ':SizeMismatch'], 'f and H must have the same length.');
end

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Passband", [], @(v) isnumeric(v) && (isempty(v) || (ismatrix(v) && size(v,2) == 2)));
addParameter(p, "Stopband", [], @(v) isnumeric(v) && (isempty(v) || (ismatrix(v) && size(v,2) == 2)));
addParameter(p, "PassbandGain", 1, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
addParameter(p, "Eps", 1e-20, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
parse(p, varargin{:});

pb = p.Results.Passband;
sb = p.Results.Stopband;
gRef = p.Results.PassbandGain;
eps0 = p.Results.Eps;

mag = abs(H);
spec = struct();

if ~isempty(pb)
    maskPb = false(size(f));
    for k = 1:size(pb,1)
        f1 = pb(k,1); f2 = pb(k,2);
        if ~(f1 <= f2)
            error([mfilename ':InvalidPassband'], 'Passband rows must satisfy f1 <= f2.');
        end
        maskPb = maskPb | (f >= f1 & f <= f2);
    end
    if ~any(maskPb)
        error([mfilename ':EmptyPassbandMask'], 'Passband mask is empty for the given frequency axis.');
    end
    dbPb = 20*log10(mag(maskPb) / gRef + eps0);
    spec.passbandMaxDb = max(dbPb);
    spec.passbandMinDb = min(dbPb);
    spec.passbandRippleDb = spec.passbandMaxDb - spec.passbandMinDb;
end

if ~isempty(sb)
    maskSb = false(size(f));
    for k = 1:size(sb,1)
        f1 = sb(k,1); f2 = sb(k,2);
        if ~(f1 <= f2)
            error([mfilename ':InvalidStopband'], 'Stopband rows must satisfy f1 <= f2.');
        end
        maskSb = maskSb | (f >= f1 & f <= f2);
    end
    if ~any(maskSb)
        error([mfilename ':EmptyStopbandMask'], 'Stopband mask is empty for the given frequency axis.');
    end
    dbSb = 20*log10(mag(maskSb) + eps0);
    spec.stopbandMaxDb = max(dbSb);      % closest to 0 dB
    spec.stopbandAttenDb = -spec.stopbandMaxDb; % positive number
end
end





