function w = sp_window_rect(N)
%SP_WINDOW_RECT Rectangular window (all ones), as a column vector.
%
%   w = sp_window_rect(N) returns an N-by-1 vector of ones.
%
%   Notes:
%   - A finite-length observation is equivalent to multiplying the signal by
%     this rectangular window.

validateattributes(N, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 1);
w = ones(N, 1);
end





