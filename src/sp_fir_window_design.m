function [h, info] = sp_fir_window_design(fs, type, fc, N, varargin)
%SP_FIR_WINDOW_DESIGN FIR design by the window method (no toolboxes).
%
%   [h,info] = sp_fir_window_design(fs,type,fc,N,...) designs an N-tap linear-phase
%   FIR filter by multiplying an ideal impulse response by a window.
%
%   Supported types:
%     - "lpf" : low-pass,  fc is scalar (Hz)
%     - "hpf" : high-pass, fc is scalar (Hz)
%     - "bpf" : band-pass, fc is [f1 f2] (Hz)
%     - "bsf" : band-stop, fc is [f1 f2] (Hz)
%
%   Name-value options:
%     "Window" (default "hann") - "rect"|"hann"|"hamming"| or an N×1 vector
%     "Scale"  (default true)   - Normalize gain (DC / Nyquist / band center)
%
%   Outputs:
%     h    - FIR coefficients (N×1, real)
%     info - struct with design details
%
%   See also: sp_fir_ideal_lpf, sp_freq_response

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 1);
type = lower(string(type));
validateattributes(N, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 4);

if mod(N, 2) ~= 1
    error([mfilename ':OddLengthRequired'], 'For simplicity, this book uses odd N (Type-I linear phase). Choose odd N.');
end

if isscalar(fc)
    validateattributes(fc, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fc', 3);
else
    validateattributes(fc, {'numeric'}, {'vector', 'real', 'finite', 'positive'}, mfilename, 'fc', 3);
    if numel(fc) ~= 2
        error([mfilename ':InvalidFc'], 'For band filters, fc must be a 2-element vector [f1 f2].');
    end
    fc = fc(:).';
end

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Window", "hann", @(v) (isstring(v) || ischar(v)) || (isnumeric(v) && isvector(v)));
addParameter(p, "Scale", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
parse(p, varargin{:});

winIn = p.Results.Window;
doScale = logical(p.Results.Scale);

% Window vector
if isstring(winIn) || ischar(winIn)
    wName = lower(string(winIn));
    switch wName
        case "rect"
            w = sp_window_rect(N);
        case "hann"
            w = sp_window_hann(N);
        case "hamming"
            w = sp_window_hamming(N);
        otherwise
            error([mfilename ':InvalidWindow'], 'Window must be "rect", "hann", "hamming", or an N×1 vector.');
    end
else
    w = winIn(:);
    wName = "custom";
    if numel(w) ~= N
        error([mfilename ':WindowSizeMismatch'], 'Custom window must have length N.');
    end
end

alpha = (N-1)/2;
d = zeros(N, 1);
d(alpha + 1) = 1; % delta[n-alpha]

switch type
    case "lpf"
        hIdeal = sp_fir_ideal_lpf(fs, fc, N);
    case "hpf"
        hIdeal = d - sp_fir_ideal_lpf(fs, fc, N);
    case "bpf"
        f1 = fc(1); f2 = fc(2);
        if ~(f1 < f2)
            error([mfilename ':InvalidBand'], 'Band edges must satisfy f1 < f2.');
        end
        hIdeal = sp_fir_ideal_lpf(fs, f2, N) - sp_fir_ideal_lpf(fs, f1, N);
    case "bsf"
        f1 = fc(1); f2 = fc(2);
        if ~(f1 < f2)
            error([mfilename ':InvalidBand'], 'Band edges must satisfy f1 < f2.');
        end
        hIdeal = d - (sp_fir_ideal_lpf(fs, f2, N) - sp_fir_ideal_lpf(fs, f1, N));
    otherwise
        error([mfilename ':InvalidType'], 'type must be "lpf","hpf","bpf","bsf".');
end

h = hIdeal .* w;

scale = 1;
if doScale
    n = (0:N-1).';
    switch type
        case "lpf"
            % DC gain H(0) = sum h[n]
            g = sum(h);
            scale = 1 / g;
        case "hpf"
            % Nyquist gain H(pi) = sum h[n] (-1)^n
            g = sum(h .* (-1).^n);
            scale = 1 / g;
        case "bpf"
            f0 = mean(fc);
            Omega0 = 2*pi*f0/fs;
            g = sum(h .* exp(-1j*Omega0*n));
            scale = 1 / abs(g);
        case "bsf"
            g = sum(h);
            scale = 1 / g;
    end
    h = real(scale * h); % keep real (linear-phase)
end

info = struct();
info.fs = fs;
info.type = type;
info.fc = fc;
info.N = N;
info.window = w;
info.windowName = wName;
info.groupDelaySamples = (N-1)/2;
info.scale = scale;
end





