function [x, t, n] = sp_complex_tone(fs, N, f0, varargin)
%SP_COMPLEX_TONE Complex exponential (complex sinusoid) signal generator.
%
%   x = sp_complex_tone(fs,N,f0) generates N samples of a complex exponential
%   at frequency f0 [Hz] sampled at fs [Hz]:
%       x[n] = A * exp(1j * (2*pi*f0*t[n] + phi)),  t[n]=n/fs.
%
%   [x,t,n] = ... also returns the time axis t [s] and the sample index n
%   (both N-by-1 column vectors).
%
%   Name-value options:
%     "Amplitude" (default 1) - Signal amplitude (real scalar)
%     "Phase"     (default 0) - Initial phase [rad] (real scalar)

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 1);
validateattributes(N, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 2);
validateattributes(f0, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, 'f0', 3);

p = inputParser;
p.FunctionName = mfilename;
p.KeepUnmatched = false;
addParameter(p, "Amplitude", 1, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite'}));
addParameter(p, "Phase", 0, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite'}));
parse(p, varargin{:});

amplitude = p.Results.Amplitude;
phase = p.Results.Phase;

[t, n] = sp_time_axis(fs, N);
x = amplitude * exp(1j * (2*pi*f0*t + phase));
end

