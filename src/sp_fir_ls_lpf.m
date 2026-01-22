function [h, info] = sp_fir_ls_lpf(fs, Fp, Fs, N, varargin)
%SP_FIR_LS_LPF Least-squares Type-I linear-phase low-pass FIR design (no toolboxes).
%
%   [h,info] = sp_fir_ls_lpf(fs,Fp,Fs,N,...) designs an odd-length (Type-I)
%   linear-phase FIR low-pass filter by weighted least squares on a frequency grid.
%
%   The design matches the real even amplitude A(Ω) in:
%     H(e^{jΩ}) = e^{-jΩα} A(Ω),   α = (N-1)/2
%   where A(Ω) is represented as a cosine series:
%     A(Ω) = c0 + 2*sum_{m=1..M} c_m cos(Ω m),   M=(N-1)/2
%
%   The coefficients h[n] are then formed by symmetry:
%     h[α]     = c0
%     h[α±m]   = c_m
%
%   Name-value options:
%     "GridN"  (default 4096) - Number of frequency samples on 0..fs/2
%     "Wp"     (default 1)    - Passband weight
%     "Ws"     (default 30)   - Stopband weight
%     "Scale"  (default true) - Normalize DC gain to 1
%
%   Notes:
%   - The transition band (Fp..Fs) is ignored in the fit (weight 0).
%   - This is a learning-oriented minimal LS design (not an equiripple design).
%
%   Inputs:
%     fs - Sampling frequency [Hz]
%     Fp - Passband edge [Hz]
%     Fs - Stopband edge [Hz]
%     N  - Number of taps (odd integer)
%
%   Outputs:
%     h    - FIR coefficients (N×1, real, symmetric)
%     info - struct with design details
%
%   See also: sp_fir_window_design, sp_filter_specs_eval, sp_freq_response

validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fs', 1);
validateattributes(Fp, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'Fp', 2);
validateattributes(Fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'Fs', 3);
validateattributes(N,  {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'N', 4);

if ~(Fp < Fs)
    error([mfilename ':InvalidEdges'], 'Edges must satisfy Fp < Fs.');
end
if Fs >= fs/2
    error([mfilename ':InvalidFs'], 'Fs must be < fs/2.');
end
if mod(N, 2) ~= 1
    error([mfilename ':OddLengthRequired'], 'This function designs Type-I linear-phase FIR: choose odd N.');
end

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, "GridN", 4096, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}));
addParameter(p, "Wp", 1, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
addParameter(p, "Ws", 30, @(v) validateattributes(v, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
addParameter(p, "Scale", true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
parse(p, varargin{:});

gridN = double(p.Results.GridN);
Wp = double(p.Results.Wp);
Ws = double(p.Results.Ws);
doScale = logical(p.Results.Scale);

fGrid = linspace(0, fs/2, gridN).';
Omega = 2*pi*fGrid/fs; % rad/sample

maskPb = (fGrid <= Fp);
maskSb = (fGrid >= Fs);

fFit = fGrid(maskPb | maskSb);
OmegaFit = Omega(maskPb | maskSb);

d = double(maskPb(maskPb | maskSb)); % 1 in PB, 0 in SB
w = ones(size(d));
w(d == 1) = Wp;
w(d == 0) = Ws;

% Build cosine-series matrix for A(Ω)
M = (N-1)/2;
mm = (0:M);
A = cos(OmegaFit * mm);   % cos(Ω m), includes m=0 column
A(:, 2:end) = 2 * A(:, 2:end);

Aw = A .* w; % row scaling
dw = d .* w;

c = Aw \ dw; % least squares

% Form symmetric h from c
alpha = M;
h = zeros(N, 1);
h(alpha + 1) = c(1);
for m = 1:M
    h(alpha + 1 + m) = c(m + 1);
    h(alpha + 1 - m) = c(m + 1);
end

scale = 1;
if doScale
    g0 = sum(h); % DC gain
    scale = 1 / g0;
    h = scale * h;
end

% Ensure real (numerical safety)
h = real(h);

info = struct();
info.fs = fs;
info.Fp = Fp;
info.Fs = Fs;
info.N = N;
info.type = "lpf";
info.method = "ls";
info.gridN = gridN;
info.weights = struct("Wp", Wp, "Ws", Ws);
info.groupDelaySamples = (N-1)/2;
info.scale = scale;
info.fFit = fFit;
end





