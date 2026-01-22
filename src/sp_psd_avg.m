function [f, Pxx, info] = sp_psd_avg(x, fs, varargin)
%SP_PSD_AVG Welch-like averaged PSD estimate (no toolboxes).
%
%   [f,Pxx] = sp_psd_avg(x,fs) estimates the power spectral density (PSD) of
%   x sampled at fs [Hz] by:
%     1) framing (overlap allowed),
%     2) windowing,
%     3) FFT,
%     4) magnitude-squared,
%     5) averaging across frames.
%
%   Name-value options:
%     "FrameLength" (default min(1024,length(x))) - Frame length [samples]
%     "Hop"         (default floor(FrameLength/2)) - Hop size [samples]
%     "Window"      (default "hann") - "rect","hann","hamming", or a vector
%     "Nfft"        (default FrameLength) - FFT length (>= FrameLength)
%     "Range"       (default "half") - "half": 0..fs/2 (real signals only)
%                                      "whole": -fs/2..fs/2 (two-sided)
%     "PadEnd"      (default false) - Pad the end so the last partial frame is included
%
%   Outputs:
%     f    - Frequency axis [Hz] (column vector)
%     Pxx  - PSD estimate [unit^2/Hz] (column vector)
%     info - Struct with details (nFrames, df, scaling, etc.)
%
%   Notes:
%   - The scaling follows the common periodogram/Welch convention:
%       P2[k] = |X[k]|^2 / (fs * sum(w.^2))
%     where X is the Nfft-point FFT of the windowed frame.
%   - When Range="half", the PSD is converted to one-sided form by doubling
%     non-DC (and non-Nyquist for even Nfft) bins.

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 2);
x = x(:);
if isempty(x)
    error([mfilename ':EmptyInput'], 'x must be non-empty.');
end

defaultFrameLength = min(1024, numel(x));
defaultFrameLength = max(8, defaultFrameLength);
defaultFrameLength = min(defaultFrameLength, numel(x));

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "FrameLength", defaultFrameLength, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}));
addParameter(p, "Hop", [], @(v) isempty(v) || (isscalar(v) && isnumeric(v) && isfinite(v) && v > 0 && mod(v,1) == 0));
addParameter(p, "Window", "hann", @(v) (isstring(v) || ischar(v)) || (isnumeric(v) && isvector(v) && ~isempty(v)));
addParameter(p, "Nfft", [], @(v) isempty(v) || (isscalar(v) && isnumeric(v) && isfinite(v) && v > 0 && mod(v,1) == 0));
addParameter(p, "Range", "half", @(v) isstring(v) || ischar(v));
addParameter(p, "PadEnd", false, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
parse(p, varargin{:});

L = p.Results.FrameLength;
hop = p.Results.Hop;
if isempty(hop)
    hop = floor(L/2);
end

nfft = p.Results.Nfft;
if isempty(nfft)
    nfft = L;
end
if nfft < L
    error([mfilename ':NfftTooSmall'], 'Nfft must be >= FrameLength.');
end

range = lower(string(p.Results.Range));
padEnd = logical(p.Results.PadEnd);

% Window
wArg = p.Results.Window;
winName = "";
if isnumeric(wArg)
    w = wArg(:);
    if numel(w) ~= L
        error([mfilename ':WindowLengthMismatch'], 'Window length must match FrameLength.');
    end
    winName = "custom";
else
    winName = lower(string(wArg));
    switch winName
        case "rect"
            w = sp_window_rect(L);
        case "hann"
            w = sp_window_hann(L);
        case "hamming"
            w = sp_window_hamming(L);
        otherwise
            error([mfilename ':InvalidWindow'], 'Window must be "rect", "hann", "hamming", or a numeric vector.');
    end
end

if range == "half" && ~isreal(x)
    error([mfilename ':HalfRangeComplex'], 'Range="half" is only supported for real x.');
end

% Framing
[X, startIdx] = sp_frame_signal(x, L, hop, "PadEnd", padEnd, "PadValue", 0);
nFrames = size(X, 2);

% Windowing (column-wise)
Xw = X .* w; % w is Lx1 (implicit expansion)

% FFT along rows (each column is a frame)
Xf = fft(Xw, nfft, 1);

winPower = sum(w.^2);
P2_frames = abs(Xf).^2 / (fs * winPower);
P2 = mean(P2_frames, 2);

df = fs / nfft;

switch range
    case "half"
        kMax = floor(nfft/2);
        Pxx = P2(1:kMax+1);
        f = (0:kMax).' * df;

        if nfft > 1
            if mod(nfft, 2) == 0
                if numel(Pxx) > 2
                    Pxx(2:end-1) = 2 * Pxx(2:end-1);
                end
            else
                if numel(Pxx) > 1
                    Pxx(2:end) = 2 * Pxx(2:end);
                end
            end
        end

    case "whole"
        Pxx = fftshift(P2);
        [f, ~] = sp_fft_axis(fs, nfft, "Range", "whole");

    otherwise
        error([mfilename ':InvalidRange'], 'Range must be "half" or "whole".');
end

if nargout >= 3
    info = struct();
    info.frameLength = L;
    info.hop = hop;
    info.nfft = nfft;
    info.df = df;
    info.nFrames = nFrames;
    info.startIdx = startIdx;
    info.window = winName;
    info.windowPower = winPower;
    info.coherentGain = mean(w);
    info.range = range;
    info.padEnd = padEnd;
end
end


