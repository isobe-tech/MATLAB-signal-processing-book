function ax = sp_plot_time(t, x, varargin)
%SP_PLOT_TIME Plot time-domain waveform with consistent labels.
%
%   ax = sp_plot_time(t,x) plots x versus t and returns the axes handle.
%
%   Name-value options:
%     "Title"     (default "")           - Plot title
%     "XLabel"    (default "Time [s]")   - X-axis label
%     "YLabel"    (default "Amplitude")  - Y-axis label
%     "NewFigure" (default true)         - Create a new figure
%     "Grid"      (default true)         - Show grid
%     "LineWidth" (default 1.2)          - Line width
%     "Part"      (default "real")       - For complex x: "real","imag","abs","phase"

t = t(:);
x = x(:);
if numel(t) ~= numel(x)
    error([mfilename ':SizeMismatch'], 't and x must have the same length.');
end

defaultPart = "real";
p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Title", "", @(v) isstring(v) || ischar(v));
addParameter(p, "XLabel", "Time [s]", @(v) isstring(v) || ischar(v));
addParameter(p, "YLabel", "Amplitude", @(v) isstring(v) || ischar(v));
addParameter(p, "NewFigure", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "Grid", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "LineWidth", 1.2, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
addParameter(p, "Part", defaultPart, @(v) isstring(v) || ischar(v));
parse(p, varargin{:});

plotTitle = string(p.Results.Title);
xLabel = string(p.Results.XLabel);
yLabel = string(p.Results.YLabel);
newFigure = logical(p.Results.NewFigure);
gridOn = logical(p.Results.Grid);
lineWidth = p.Results.LineWidth;
part = validatestring(p.Results.Part, ["real", "imag", "abs", "phase"], mfilename, "Part");

switch part
    case 'real'
        y = real(x);
    case 'imag'
        y = imag(x);
    case 'abs'
        y = abs(x);
    case 'phase'
        y = unwrap(angle(x));
end

if newFigure
    figure;
end

ax = gca;
plot(ax, t, y, "LineWidth", lineWidth);
xlabel(ax, xLabel);
ylabel(ax, yLabel);
title(ax, plotTitle);
if gridOn
    grid(ax, "on");
else
    grid(ax, "off");
end
end
