function y = sp_iir2(x, b, a, varargin)
%SP_IIR2 Second-order IIR (biquad) implemented by a difference equation (no toolboxes).
%
%   y = sp_iir2(x, b, a) filters input x by the 2nd-order IIR transfer function:
%     H(z) = (b0 + b1 z^{-1} + b2 z^{-2}) / (a0 + a1 z^{-1} + a2 z^{-2})
%
%   where:
%     b = [b0 b1 b2] and a = [a0 a1 a2] (row/column ok).
%
%   The time-domain recursion (Direct Form I) is:
%     y[n] = ( b0 x[n] + b1 x[n-1] + b2 x[n-2]
%              - a1 y[n-1] - a2 y[n-2] ) / a0
%
%   Name-value options:
%     "Structure"  (default "df1")   - "df1" (direct-form I) or "df2" (direct-form II)
%     "Arithmetic" (default "double")- "double" or "single" (internal arithmetic)
%     "XPast"      (default [0;0])   - For df1: [x[n-1]; x[n-2]]
%     "YPast"      (default [0;0])   - For df1: [y[n-1]; y[n-2]]
%     "WPast"      (default [0;0])   - For df2: [w[n-1]; w[n-2]] (canonical state)
%
%   Output:
%     y - Output signal (N×1 column vector)
%
%   Notes:
%   - This function is intentionally limited to 2nd order to keep the
%     implementation readable for learning and for "biquad / SOS" discussions.
%   - For verification/comparison, MATLAB's filter(b,a,x) can be used.
%
%   See also: sp_freq_response, sp_group_delay

validateattributes(x, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'x', 1);
validateattributes(b, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'b', 2);
validateattributes(a, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'a', 3);

b = b(:);
a = a(:);
if numel(b) ~= 3 || numel(a) ~= 3
    error([mfilename ':InvalidOrder'], 'b and a must be length-3 vectors: [b0 b1 b2], [a0 a1 a2].');
end
if a(1) == 0
    error([mfilename ':InvalidA0'], 'a(1) (=a0) must be non-zero.');
end

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Structure", "df1", @(v) isstring(v) || ischar(v));
addParameter(p, "Arithmetic", "double", @(v) isstring(v) || ischar(v));
addParameter(p, "XPast", [0; 0], @(v) isnumeric(v) && isvector(v) && numel(v) == 2);
addParameter(p, "YPast", [0; 0], @(v) isnumeric(v) && isvector(v) && numel(v) == 2);
addParameter(p, "WPast", [0; 0], @(v) isnumeric(v) && isvector(v) && numel(v) == 2);
parse(p, varargin{:});

structure = lower(string(p.Results.Structure));
arith = lower(string(p.Results.Arithmetic));

switch arith
    case "double"
        cls = "double";
    case "single"
        cls = "single";
    otherwise
        error([mfilename ':InvalidArithmetic'], 'Arithmetic must be "double" or "single".');
end

x = cast(x(:), cls);
b = cast(b, cls);
a = cast(a, cls);

b0 = b(1); b1 = b(2); b2 = b(3);
a0 = a(1); a1 = a(2); a2 = a(3);

N = numel(x);
y = zeros(N, 1, cls);

switch structure
    case "df1"
        xPast = cast(p.Results.XPast(:), cls);
        yPast = cast(p.Results.YPast(:), cls);
        x1 = xPast(1); x2 = xPast(2);
        y1 = yPast(1); y2 = yPast(2);

        for n = 1:N
            yc = (b0*x(n) + b1*x1 + b2*x2 - a1*y1 - a2*y2) / a0;
            y(n) = yc;
            x2 = x1; x1 = x(n);
            y2 = y1; y1 = yc;
        end

    case "df2"
        wPast = cast(p.Results.WPast(:), cls);
        w1 = wPast(1); w2 = wPast(2);

        for n = 1:N
            w0 = (x(n) - a1*w1 - a2*w2) / a0;
            yc = b0*w0 + b1*w1 + b2*w2;
            y(n) = yc;
            w2 = w1; w1 = w0;
        end

    otherwise
        error([mfilename ':InvalidStructure'], 'Structure must be "df1" or "df2".');
end
end





