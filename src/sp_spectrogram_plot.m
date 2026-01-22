function ax = sp_spectrogram_plot(t, f, P, varargin)
%SP_SPECTROGRAM_PLOT Plot a spectrogram-like image with consistent conventions.
%
%   ax = sp_spectrogram_plot(t,f,P) visualizes a time-frequency matrix P
%   with time axis t [s] and frequency axis f [Hz].
%
%   P is typically a PSD matrix [unit^2/Hz] (as returned by sp_stft with
%   Output="psd"). This function plots 10*log10(P) by default.
%
%   Name-value options:
%     "Title"            (default "")     - Title
%     "Db"               (default true)   - Plot in dB
%     "DbOffset"         (default 1e-20)  - Added before log to avoid -Inf
%     "DynamicRangeDb"   (default 80)     - Clip to [max-DR, max]
%     "NewFigure"        (default true)   - Create a new figure
%     "Colormap"         (default "parula") - Colormap name
%     "ColorbarLabel"    (default "PSD [dB/Hz]") - Label for colorbar
%
%   Notes:
%   - Axis is set to "xy" so that low frequencies appear at the bottom.

t = t(:);
f = f(:);
validateattributes(P, {'numeric'}, {'2d', 'nonempty', 'real', 'finite', 'nonnegative'}, mfilename, 'P', 3);

if size(P, 1) ~= numel(f) || size(P, 2) ~= numel(t)
    error([mfilename ':SizeMismatch'], 'P must be size numel(f)-by-numel(t).');
end

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Title", "", @(v) isstring(v) || ischar(v));
addParameter(p, "Db", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "DbOffset", 1e-20, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'}));
addParameter(p, "DynamicRangeDb", 80, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
addParameter(p, "NewFigure", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "Colormap", "parula", @(v) isstring(v) || ischar(v));
addParameter(p, "ColorbarLabel", "PSD [dB/Hz]", @(v) isstring(v) || ischar(v));
parse(p, varargin{:});

plotTitle = string(p.Results.Title);
useDb = logical(p.Results.Db);
dbOffset = p.Results.DbOffset;
dr = p.Results.DynamicRangeDb;
newFigure = logical(p.Results.NewFigure);
cmap = string(p.Results.Colormap);
cbLabel = string(p.Results.ColorbarLabel);

if useDb
    Z = 10 * log10(P + dbOffset);
    zMax = max(Z(:));
    Z = max(Z, zMax - dr);
else
    Z = P;
end

if newFigure
    figure;
end
ax = gca;
imagesc(ax, t, f, Z);
axis(ax, "xy");
xlabel(ax, "Time [s]");
ylabel(ax, "Frequency [Hz]");
title(ax, plotTitle);
colormap(ax, cmap);
cb = colorbar(ax);
cb.Label.String = cbLabel;
try
    cb.Color = [0 0 0];
    cb.Label.Color = [0 0 0];
catch
end
end




