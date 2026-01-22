function [t, f, S, info] = sp_stft(x, fs, varargin)
%SP_STFT Short-time Fourier transform / spectrogram building block (no toolboxes).
%
%   [t,f,S] = sp_stft(x,fs) computes the STFT of x sampled at fs [Hz] using
%   framing + windowing + FFT. The output S is either complex STFT coefficients
%   or per-frame PSD, depending on "Output".
%
%   Name-value options:
%     "FrameLength"   (default min(1024,length(x))) - Frame length L [samples]
%     "Hop"           (default floor(L/2))          - Hop size R [samples]
%     "Window"        (default "hann")              - "rect","hann","hamming", or vector
%     "Nfft"          (default L)                   - FFT length (>= L)
%     "Range"         (default "half")              - "half": 0..fs/2 (real x only)
%                                                     "whole": -fs/2..fs/2 (fftshift)
%     "PadEnd"        (default false)               - Pad the end so the last frame is included
%     "TimeReference" (default "center")            - "start" or "center" time for each frame
%     "Output"        (default "psd")               - "psd" or "complex"
%
%   Outputs:
%     t    - Frame times [s] (column vector; according to TimeReference)
%     f    - Frequency axis [Hz] (column vector; according to Range)
%     S    - STFT output:
%            Output="complex": complex STFT coefficients
%            Output="psd"    : PSD per frame [unit^2/Hz] (real, >=0)
%     info - Struct with parameters and helpful metadata (startIdx, df, etc.)
%
%   Notes:
%   - Output="psd" uses the same scaling as Chapter 7:
%       P2[k] = |X[k]|^2 / (fs * sum(w.^2))
%     and one-sided conversion doubles non-DC (and non-Nyquist if even Nfft).

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
addParameter(p, "TimeReference", "center", @(v) isstring(v) || ischar(v));
addParameter(p, "Output", "psd", @(v) isstring(v) || ischar(v));
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
timeRef = lower(string(p.Results.TimeReference));
outMode = lower(string(p.Results.Output));

if range == "half" && ~isreal(x)
    error([mfilename ':HalfRangeComplex'], 'Range="half" is only supported for real x.');
end

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

% Framing
[X, startIdx] = sp_frame_signal(x, L, hop, "PadEnd", padEnd, "PadValue", 0);
nFrames = size(X, 2);

% Time axis
switch timeRef
    case "start"
        nRef = double(startIdx) - 1;
    case "center"
        nRef = double(startIdx) - 1 + (L-1)/2;
    otherwise
        error([mfilename ':InvalidTimeReference'], 'TimeReference must be "start" or "center".');
end
t = nRef(:) / fs;

% Windowing and FFT (each column is a frame)
Xw = X .* w; % implicit expansion
Xf = fft(Xw, nfft, 1);

df = fs / nfft;

% Convert to requested output
switch outMode
    case "complex"
        Z = Xf;
    case "psd"
        winPower = sum(w.^2);
        Z = abs(Xf).^2 / (fs * winPower);
    otherwise
        error([mfilename ':InvalidOutput'], 'Output must be "complex" or "psd".');
end

% Frequency axis and range selection
switch range
    case "half"
        kMax = floor(nfft/2);
        S = Z(1:kMax+1, :);
        f = (0:kMax).' * df;

        if outMode == "psd" && nfft > 1
            if mod(nfft, 2) == 0
                if size(S, 1) > 2
                    S(2:end-1, :) = 2 * S(2:end-1, :);
                end
            else
                if size(S, 1) > 1
                    S(2:end, :) = 2 * S(2:end, :);
                end
            end
        end

    case "whole"
        S = fftshift(Z, 1);
        [f, ~] = sp_fft_axis(fs, nfft, "Range", "whole");

    otherwise
        error([mfilename ':InvalidRange'], 'Range must be "half" or "whole".');
end

if nargout >= 4
    info = struct();
    info.frameLength = L;
    info.hop = hop;
    info.nfft = nfft;
    info.df = df;
    info.nFrames = nFrames;
    info.startIdx = startIdx;
    info.window = winName;
    info.windowPower = sum(w.^2);
    info.coherentGain = mean(w);
    info.range = range;
    info.padEnd = padEnd;
    info.timeReference = timeRef;
    info.output = outMode;
end
end





