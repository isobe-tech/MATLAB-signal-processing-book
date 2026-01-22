function y = sp_clip(x, lo, hi)
%SP_CLIP Saturate values to the closed interval [lo, hi].
%
%   y = sp_clip(x, lo, hi) clips x to the range [lo, hi] elementwise.
%   This utility is used in Chapter 2 to model clipping (saturation).
%
%   Inputs:
%     x  - Signal (vector). It will be treated as a column vector.
%     lo - Lower bound (real scalar)
%     hi - Upper bound (real scalar, hi >= lo)
%
%   Output:
%     y  - Clipped signal (column vector)

validateattributes(lo, {'numeric'}, {'scalar', 'real', 'finite'}, mfilename, 'lo', 2);
validateattributes(hi, {'numeric'}, {'scalar', 'real', 'finite'}, mfilename, 'hi', 3);
if hi < lo
    error([mfilename ':InvalidRange'], 'hi must be >= lo.');
end

x = x(:);
y = min(max(x, lo), hi);
end

