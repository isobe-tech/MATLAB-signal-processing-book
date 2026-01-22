function [x, t, xClean, w, noiseStd] = sp_sine_with_noise(fs, N, f0, snrDb, varargin)
%SP_SINE_WITH_NOISE Sine wave with additive white noise (SNR specified).
%
%   x = sp_sine_with_noise(fs,N,f0,snrDb) generates N samples of a sine wave
%   at frequency f0 [Hz] sampled at fs [Hz], and adds white Gaussian noise
%   such that the SNR (dB) is snrDb, where power is defined as mean(|x|^2).
%
%   [x,t,xClean,w,noiseStd] = ... also returns the time axis t [s], the clean
%   signal xClean, the noise w, and the noise standard deviation. For complex
%   noise, noiseStd is the standard deviation per real/imag component.
%
%   Name-value options:
%     "Amplitude" (default 1)    - Signal amplitude (real scalar)
%     "Phase"     (default 0)    - Initial phase [rad] (real scalar)
%     "Complex"   (default false)- true -> complex exponential
%     "Seed"      (default [])   - RNG seed for reproducibility (integer)

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 1);
validateattributes(N, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 2);
validateattributes(f0, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, 'f0', 3);
validateattributes(snrDb, {'numeric'}, {'scalar', 'real'}, mfilename, 'snrDb', 4);
if isnan(snrDb) || (isinf(snrDb) && snrDb < 0)
    error([mfilename ':snrDb'], 'snrDb must be a real scalar (finite or +Inf).');
end

p = inputParser;
p.FunctionName = mfilename;
p.KeepUnmatched = false;
addParameter(p, "Amplitude", 1, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite'}));
addParameter(p, "Phase", 0, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite'}));
addParameter(p, "Complex", false, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "Seed", [], @(v) isempty(v) || (isscalar(v) && isnumeric(v) && isfinite(v) && v >= 0 && mod(v, 1) == 0));
parse(p, varargin{:});

amplitude = p.Results.Amplitude;
phase = p.Results.Phase;
isComplex = logical(p.Results.Complex);
seed = p.Results.Seed;

if ~isempty(seed)
    oldState = rng;
    cleanupRng = onCleanup(@() rng(oldState)); %#ok<NASGU>
    rng(seed, "twister");
end

[t, ~] = sp_time_axis(fs, N);

if isComplex
    xClean = amplitude * exp(1j * (2*pi*f0*t + phase));
else
    xClean = amplitude * sin(2*pi*f0*t + phase);
end

if isinf(snrDb)
    w = zeros(N, 1);
    noiseStd = 0;
    x = xClean;
    return;
end

signalPower = mean(abs(xClean).^2);
noisePower = signalPower / (10^(snrDb/10));

if isComplex
    noiseStd = sqrt(noisePower/2);
    w = noiseStd * (randn(N, 1) + 1j*randn(N, 1));
else
    noiseStd = sqrt(noisePower);
    w = noiseStd * randn(N, 1);
end

x = xClean + w;
end
