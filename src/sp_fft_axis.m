function [f, k] = sp_fft_axis(fs, Nfft, varargin)
%SP_FFT_AXIS Frequency axis for FFT/DFT bins (no toolboxes).
%
%   [f,k] = sp_fft_axis(fs,Nfft) returns a frequency axis f [Hz] for DFT bins
%   with sampling frequency fs and transform length Nfft.
%
%   Name-value options:
%     "Range" (default "half")
%       - "half": 0 .. fs/2  (one-sided, k = 0..floor(Nfft/2))
%       - "whole": -fs/2 .. fs/2 (two-sided, centered at 0; for fftshift)
%       - "raw": 0 .. fs*(Nfft-1)/Nfft (unshifted full bins)
%
%   Outputs:
%     f - Frequency axis [Hz] (column vector)
%     k - Bin indices corresponding to f (column vector)

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 1);
validateattributes(Nfft, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'Nfft', 2);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "Range", "half", @(v) isstring(v) || ischar(v));
parse(p, varargin{:});

range = lower(string(p.Results.Range));
switch range
    case "half"
        k = (0:floor(Nfft/2)).';
        f = k * (fs / Nfft);
    case "whole"
        if mod(Nfft, 2) == 0
            k = (-Nfft/2:Nfft/2-1).';
        else
            k = (-(Nfft-1)/2:(Nfft-1)/2).';
        end
        f = k * (fs / Nfft);
    case "raw"
        k = (0:Nfft-1).';
        f = k * (fs / Nfft);
    otherwise
        error([mfilename ':InvalidRange'], 'Range must be "half", "whole", or "raw".');
end
end

