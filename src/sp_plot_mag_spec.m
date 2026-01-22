function ax = sp_plot_mag_spec(fs, x, varargin)
%SP_PLOT_MAG_SPEC Quick-look magnitude spectrum plot (FFT-based).
%
%   ax = sp_plot_mag_spec(fs,x) plots the magnitude spectrum of x sampled at
%   fs [Hz]. This function is intended as a quick verification tool; detailed
%   discussions (windowing, leakage, amplitude calibration) are handled in
%   later chapters.
%
%   Name-value options:
%     "Title"    (default "")     - Plot title
%     "Nfft"     (default length(x)) - FFT length (>= length(x))
%     "Db"       (default true)   - Plot in dB (20*log10)
%     "OneSided" (default isreal(x)) - One-sided spectrum for real signals
%     "NewFigure"(default true)   - Create a new figure
%     "Grid"     (default true)   - Show grid

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 1);
x = x(:);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Title", "", @(v) isstring(v) || ischar(v));
addParameter(p, "Nfft", numel(x), @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}));
addParameter(p, "Db", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "OneSided", isreal(x), @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "NewFigure", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "Grid", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
parse(p, varargin{:});

plotTitle = string(p.Results.Title);
nfft = p.Results.Nfft;
useDb = logical(p.Results.Db);
oneSided = logical(p.Results.OneSided);
newFigure = logical(p.Results.NewFigure);
gridOn = logical(p.Results.Grid);

N = numel(x);
if nfft < N
    error([mfilename ':NfftTooSmall'], 'Nfft must be >= length(x).');
end
if oneSided && ~isreal(x)
    error([mfilename ':OneSidedComplex'], 'OneSided spectrum is only supported for real x.');
end

X = fft(x, nfft);
mag = abs(X) / N;

if oneSided
    kMax = floor(nfft/2);
    mag = mag(1:kMax+1);
    f = (0:kMax).' * (fs / nfft);

    if nfft > 1
        if mod(nfft, 2) == 0
            if numel(mag) > 2
                mag(2:end-1) = 2 * mag(2:end-1);
            end
        else
            if numel(mag) > 1
                mag(2:end) = 2 * mag(2:end);
            end
        end
    end
else
    X = fftshift(X);
    mag = abs(X) / N;
    if mod(nfft, 2) == 0
        k = (-nfft/2:nfft/2-1).';
    else
        k = (-(nfft-1)/2:(nfft-1)/2).';
    end
    f = k * (fs / nfft);
end

if useDb
    y = 20 * log10(mag + eps);
    yLabel = "Magnitude [dB]";
else
    y = mag;
    yLabel = "Magnitude";
end

if newFigure
    figure;
end
ax = gca;
plot(ax, f, y, "LineWidth", 1.2);
xlabel(ax, "Frequency [Hz]");
ylabel(ax, yLabel);
title(ax, plotTitle);
if gridOn
    grid(ax, "on");
else
    grid(ax, "off");
end
end
