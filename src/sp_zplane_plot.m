function ax = sp_zplane_plot(b, a, varargin)
%SP_ZPLANE_PLOT Plot zeros and poles on the z-plane with the unit circle (no toolboxes).
%
%   ax = sp_zplane_plot(b,a) plots zeros (o) and poles (x) of
%     H(z) = B(z)/A(z)
%   where coefficient vectors b and a represent polynomials in z^{-1}:
%     B(z) = b0 + b1 z^{-1} + ...,
%     A(z) = a0 + a1 z^{-1} + ...
%
%   For plotting, the pole/zero locations are the roots of the corresponding
%   polynomials in z (the same coefficient vector can be used as a polynomial
%   in z after multiplying by a power of z).
%
%   Name-value options:
%     "Title"        (default "")     - Title
%     "NewFigure"    (default true)   - Create a new figure
%     "Grid"         (default true)   - Show grid
%     "AxisLimit"    (default 1.6)    - Axis limit for both real/imag
%     "MarkerSize"   (default 7)      - Marker size
%     "LineWidth"    (default 1.2)    - Line width
%     "ShowUnitCircle" (default true) - Draw the unit circle
%
%   Output:
%     ax - Axes handle
%
%   See also: roots

validateattributes(b, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'b', 1);
validateattributes(a, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'a', 2);

b = b(:);
a = a(:);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Title", "", @(v) isstring(v) || ischar(v));
addParameter(p, "NewFigure", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "Grid", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "AxisLimit", 1.6, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
addParameter(p, "MarkerSize", 7, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
addParameter(p, "LineWidth", 1.2, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
addParameter(p, "ShowUnitCircle", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "LegendLocation", "northwest", @(v) isstring(v) || ischar(v));
parse(p, varargin{:});

plotTitle = string(p.Results.Title);
newFigure = logical(p.Results.NewFigure);
gridOn = logical(p.Results.Grid);
axisLim = p.Results.AxisLimit;
ms = p.Results.MarkerSize;
lw = p.Results.LineWidth;
showUc = logical(p.Results.ShowUnitCircle);
legendLoc = string(p.Results.LegendLocation);

if newFigure
    figure;
end
ax = gca;
hold(ax, "on");

% Unit circle
if showUc
    th = linspace(0, 2*pi, 361);
    plot(ax, cos(th), sin(th), "k-", "LineWidth", 1.0);
end

% Axes lines
plot(ax, [-axisLim axisLim], [0 0], "k:", "LineWidth", 1.0);
plot(ax, [0 0], [-axisLim axisLim], "k:", "LineWidth", 1.0);

% Zeros and poles
if numel(b) <= 1
    z = [];
else
    z = roots(b);
end
if numel(a) <= 1
    pz = [];
else
    pz = roots(a);
end

if ~isempty(z)
    plot(ax, real(z), imag(z), "o", "MarkerSize", ms, "LineWidth", lw);
end
if ~isempty(pz)
    plot(ax, real(pz), imag(pz), "x", "MarkerSize", ms, "LineWidth", lw);
end

xlabel(ax, "Re\{z\}");
ylabel(ax, "Im\{z\}");
title(ax, plotTitle);
axis(ax, "equal");
xlim(ax, [-axisLim axisLim]);
ylim(ax, [-axisLim axisLim]);

if gridOn
    grid(ax, "on");
else
    grid(ax, "off");
end

legendEntries = {};
if showUc
    legendEntries{end+1} = "unit circle";
end
legendEntries{end+1} = "axes";
if ~isempty(z)
    legendEntries{end+1} = "zeros";
end
if ~isempty(pz)
    legendEntries{end+1} = "poles";
end
% Keep legend compact: only show when something exists besides guides.
if (~isempty(z) || ~isempty(pz)) && numel(legendEntries) >= 3
    % Place legend inside by default to avoid shrinking the axes area in tiledlayout.
    legend(ax, legendEntries, "Location", legendLoc);
end
hold(ax, "off");
end




