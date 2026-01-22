function h = sp_fir_ideal_lpf(fs, fc, N)
%SP_FIR_IDEAL_LPF Ideal (infinite-length) LPF impulse response sampled on a finite grid.
%
%   h = sp_fir_ideal_lpf(fs,fc,N) returns an N-by-1 vector h that corresponds to
%   the ideal low-pass frequency response:
%     H_d(e^{jΩ}) = 1  (|Ω| <= Ωc), 0 otherwise  (2π-periodic)
%   sampled around the center index (N-1)/2.
%
%   The closed-form impulse response is:
%     h[n] = sin(Ωc (n-α)) / (π (n-α)),   n ≠ α
%     h[α] = Ωc / π
%   where α = (N-1)/2 and Ωc = 2π fc / fs.
%
%   Notes:
%   - This ideal impulse response is theoretically infinite. The vector h is
%     the "centered slice" used as the starting point for the window method.
%   - For a strict Type-I linear-phase FIR (odd length), choose odd N so that
%     α is an integer sample.
%
%   Inputs:
%     fs - Sampling frequency [Hz]
%     fc - Cutoff frequency [Hz] (0 < fc < fs/2)
%     N  - Number of taps (positive integer)
%
%   Output:
%     h  - N×1 column vector (real)

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 1);
validateattributes(fc, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fc', 2);
validateattributes(N,  {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 3);

if fc >= fs/2
    error([mfilename ':InvalidFc'], 'fc must satisfy 0 < fc < fs/2.');
end

OmegaC = 2*pi*fc/fs; % rad/sample

n = (0:N-1).';
alpha = (N-1)/2;
m = n - alpha;

h = zeros(N, 1);
idx0 = (m == 0);
h(idx0) = OmegaC / pi;
h(~idx0) = sin(OmegaC * m(~idx0)) ./ (pi * m(~idx0));
end





