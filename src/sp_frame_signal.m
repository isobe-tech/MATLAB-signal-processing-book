function [X, startIdx] = sp_frame_signal(x, frameLength, hop, varargin)
%SP_FRAME_SIGNAL Slice a signal into overlapping frames (columns).
%
%   X = sp_frame_signal(x, frameLength, hop) returns a matrix X of size
%   frameLength-by-nFrames, where each column is a frame extracted from x.
%   The first frame starts at index 1 (MATLAB indexing) and subsequent frames
%   start hop samples later.
%
%   [X, startIdx] = sp_frame_signal(...) also returns the starting indices
%   (1-based) of each frame as an nFrames-by-1 vector.
%
%   Name-value options:
%     "PadEnd"   (default false) - If true, include the last partial frame by
%                                 padding beyond the end of x.
%     "PadValue" (default 0)     - Padding value used when PadEnd is true.
%
%   Notes:
%   - x is treated as a column vector.
%   - This function is intended as a small building block for Welch PSD and
%     STFT implementations without toolboxes.

x = x(:);
if isempty(x)
    error([mfilename ':EmptyInput'], 'x must be non-empty.');
end

validateattributes(frameLength, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'frameLength', 2);
validateattributes(hop, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'hop', 3);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "PadEnd", false, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
addParameter(p, "PadValue", 0, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite'}));
parse(p, varargin{:});

padEnd = logical(p.Results.PadEnd);
padValue = p.Results.PadValue;

N = numel(x);
L = frameLength;

if ~padEnd
    if N < L
        error([mfilename ':TooShort'], 'length(x) must be >= frameLength when PadEnd is false.');
    end
    nFrames = floor((N - L) / hop) + 1;
else
    % Include last partial frame by padding.
    nFrames = max(1, ceil((N - L) / hop) + 1);
end

startIdx = (0:nFrames-1).' * hop + 1;
idx = (0:L-1).' + startIdx.'; % L-by-nFrames

if ~padEnd
    X = x(idx);
else
    padVal = cast(padValue, 'like', x);
    X = padVal + zeros(L, nFrames, 'like', x);
    mask = idx <= N;
    X(mask) = x(idx(mask));
end
end





