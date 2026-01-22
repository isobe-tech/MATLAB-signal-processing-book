function sp_book_style(varargin)
%SP_BOOK_STYLE Set MATLAB graphics defaults for this book (light theme).
%
% This book embeds figures into LaTeX, so we prefer:
% - White figure/axes backgrounds
% - Black text/ticks
% - Readable font sizes (for A4, width=\linewidth figures)
%
% Usage:
%   sp_book_style();                       % defaults
%   sp_book_style("FontSize", 14);         % override
%   sp_book_style("LineWidth", 1.4);       % override
%
% Notes:
% - This function sets *root defaults* (groot), so call it near the top of
%   each chapter script before creating figures.
%

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "FontSize", 14, @(v) validateattributes(v, {'numeric'}, {'scalar','real','finite','positive'}));
addParameter(p, "LineWidth", 1.2, @(v) validateattributes(v, {'numeric'}, {'scalar','real','finite','positive'}));
addParameter(p, "AxesLineWidth", 0.75, @(v) validateattributes(v, {'numeric'}, {'scalar','real','finite','positive'}));
parse(p, varargin{:});

fs = double(p.Results.FontSize);
lw = double(p.Results.LineWidth);
alw = double(p.Results.AxesLineWidth);

% Colors
white = [1 1 1];
black = [0 0 0];
gridc = [0.85 0.85 0.85];
minorGridc = [0.90 0.90 0.90];

set(groot, ...
    "defaultFigureColor", white, ...
    "defaultAxesColor", white, ...
    "defaultAxesXColor", black, ...
    "defaultAxesYColor", black, ...
    "defaultAxesZColor", black, ...
    "defaultTextColor", black, ...
    "defaultAxesGridColor", gridc, ...
    "defaultAxesMinorGridColor", minorGridc, ...
    "defaultAxesFontSize", fs, ...
    "defaultTextFontSize", fs, ...
    "defaultLegendColor", white, ...
    "defaultLegendTextColor", black, ...
    "defaultLegendEdgeColor", black, ...
    "defaultLegendFontSize", fs, ...
    "defaultLineLineWidth", lw, ...
    "defaultAxesLineWidth", alw ...
);

% Force TeX interpreter so backslash sequences like "\alpha" render reliably
% even if a user's MATLAB preferences changed the defaults.
try
    set(groot, "defaultTextInterpreter", "tex");
catch
end
try
    set(groot, "defaultLegendInterpreter", "tex");
catch
end
try
    set(groot, "defaultAxesTickLabelInterpreter", "tex");
catch
end
try
    set(groot, "defaultColorbarTickLabelInterpreter", "tex");
catch
end

% Ensure export/print keeps the chosen colors (avoid unexpected inversion).
try
    set(groot, "defaultFigureInvertHardcopy", "off");
catch
    % Older MATLAB versions may not have this property; ignore.
end

% Avoid exporting the axes toolbar overlay (introduced in newer MATLAB versions).
try
    set(groot, "defaultAxesToolbarVisible", "off");
catch
end

end
