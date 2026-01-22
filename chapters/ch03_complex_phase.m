%% 第3章 複素指数と位相：周波数表現の基礎
% このスクリプトは「Run All」で図まで再現できる形にする。

clear; close all; clc;

% ルート（bookRoot）を推定して src/ を path に追加
thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));
sp_book_style(); % White background + readable fonts for LaTeX

% 出力先（図）
figDir = fullfile(bookRoot, "figs");
if ~exist(figDir, "dir")
    mkdir(figDir);
end

%% Parameters
fs = 8000;     % [Hz]
N = 2048;      % samples
f0 = 440;      % [Hz]
phi0 = 0.3*pi; % [rad]

omega = 2*pi*f0/fs; % [rad/sample]
fprintf("omega = 2*pi*f0/fs = %.6g rad/sample\n", omega);

%% 1) Complex exponential as rotation: real/imag vs complex plane
[x, t] = sp_complex_tone(fs, N, f0, "Phase", phi0, "Amplitude", 1);

idx = (t <= 0.01); % show first 10 ms
idxC = (1:260).';  % complex-plane segment

figure;
tiledlayout(2, 2, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(t(idx), real(x(idx)), "-", "LineWidth", 1.2);
grid on;
xlabel("Time [s]");
ylabel("Re{x}");
title("Real part");

nexttile;
plot(t(idx), imag(x(idx)), "-", "LineWidth", 1.2);
grid on;
xlabel("Time [s]");
ylabel("Im{x}");
title("Imag part");

nexttile([1 2]);
plot(real(x(idxC)), imag(x(idxC)), "o-", "LineWidth", 1.0, "MarkerSize", 2.5);
axis equal;
grid on;
xlabel("Re");
ylabel("Im");
title("Complex plane trajectory (rotation)");

exportgraphics(gcf, fullfile(figDir, "ch03_complex_rotation.png"), "Resolution", 200, "BackgroundColor", "white");

%% 2) Phase shift: waveform shift + time-shift equivalence (complex exponential)
phaseList = [0, pi/2, pi];
figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
hold on;
for k = 1:numel(phaseList)
    phi = phaseList(k);
    xk = sp_complex_tone(fs, N, f0, "Phase", phi);
    plot(t(idx), real(xk(idx)), "LineWidth", 1.1);
end
grid on;
xlabel("Time [s]");
ylabel("Re{x}");
legend("\phi=0", "\phi=\pi/2", "\phi=\pi", "Location", "southwest");
title("Changing phase moves the waveform");

% time shift equivalence on the overlapping index range
n0 = 80; % samples
x0 = sp_complex_tone(fs, N, f0, "Phase", 0);
lhs = x0(1:end-n0);                          % x[n - n0] (finite segment)
rhs = exp(-1j*omega*n0) * x0((n0+1):end);    % e^{-j*omega*n0} x[n]
relErr = norm(lhs - rhs) / max(norm(lhs), 1);
fprintf("time-shift equivalence (n0=%d): relErr=%.3g\n", n0, relErr);
assert(relErr < 1e-12, "Time-shift equivalence check failed.");

nexttile;
plot(real(lhs(1:200)), imag(lhs(1:200)), "o-", "LineWidth", 1.0, "MarkerSize", 2.2);
hold on;
plot(real(rhs(1:200)), imag(rhs(1:200)), "x-", "LineWidth", 1.0, "MarkerSize", 2.2);
axis equal;
grid on;
xlabel("Re");
ylabel("Im");
legend("x[n-n0]", "e^{-j\Omega n0} x[n]", "Location", "best");
title("Time shift equals constant phase rotation (complex exponential)");

exportgraphics(gcf, fullfile(figDir, "ch03_phase_shift.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3) Negative frequency: rotation direction and conjugate relation
xPlus = sp_complex_tone(fs, N, f0, "Phase", 0);
xMinus = conj(xPlus); % exp(-j*...) since t is real

figure;
tiledlayout(2, 2, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(t(idx), real(xPlus(idx)), "-", "LineWidth", 1.2);
hold on;
plot(t(idx), real(xMinus(idx)), "--", "LineWidth", 1.2);
grid on;
xlabel("Time [s]");
ylabel("Re{x}");
legend("+f", "-f", "Location", "southwest");
title("Real parts are identical");

nexttile;
plot(t(idx), imag(xPlus(idx)), "-", "LineWidth", 1.2);
hold on;
plot(t(idx), imag(xMinus(idx)), "--", "LineWidth", 1.2);
grid on;
xlabel("Time [s]");
ylabel("Im{x}");
legend("+f", "-f", "Location", "southwest");
title("Imag parts have opposite sign");

nexttile;
plot(real(xPlus(idxC)), imag(xPlus(idxC)), "o-", "LineWidth", 1.0, "MarkerSize", 2.5);
axis equal;
grid on;
xlabel("Re");
ylabel("Im");
title("+f : counterclockwise");

nexttile;
plot(real(xMinus(idxC)), imag(xMinus(idxC)), "o-", "LineWidth", 1.0, "MarkerSize", 2.5);
axis equal;
grid on;
xlabel("Re");
ylabel("Im");
title("-f : clockwise");

exportgraphics(gcf, fullfile(figDir, "ch03_negfreq.png"), "Resolution", 200, "BackgroundColor", "white");

%% 3.5) Real sine vs complex exponential: why real signals have ±f pairs
xC = sp_complex_tone(fs, N, f0, "Phase", 0);
xR = real(xC);

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
sp_plot_mag_spec(fs, xC, "NewFigure", false, "Db", true, "OneSided", false, ...
    "Title", "Complex exponential: one-sided component (two-sided plot)");
xlim([-2000, 2000]);

nexttile;
sp_plot_mag_spec(fs, xR, "NewFigure", false, "Db", true, "OneSided", false, ...
    "Title", "Real cosine: ±f pair (two-sided plot)");
xlim([-2000, 2000]);

exportgraphics(gcf, fullfile(figDir, "ch03_real_vs_complex_spectrum.png"), "Resolution", 200, "BackgroundColor", "white");

%% 4) Two tones: predict "two peaks" and check by quick-look spectrum
f1 = 440;
f2 = 880;
[x1, t] = sp_complex_tone(fs, N, f1, "Phase", 0.1*pi, "Amplitude", 1); %#ok<ASGLU>
x2 = sp_complex_tone(fs, N, f2, "Phase", -0.4*pi, "Amplitude", 0.6);
xSum = x1 + x2;

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(t(idx), real(xSum(idx)), "-", "LineWidth", 1.2);
grid on;
xlabel("Time [s]");
ylabel("Re{x}");
title(sprintf("Two tones in time domain: f1=%g Hz, f2=%g Hz", f1, f2));

nexttile;
sp_plot_mag_spec(fs, xSum, "NewFigure", false, "Db", true, "Title", "Magnitude spectrum (two-sided, quick look)");

exportgraphics(gcf, fullfile(figDir, "ch03_two_tones.png"), "Resolution", 200, "BackgroundColor", "white");

%% 5) Orthogonality: DFT-like complex exponentials (inner products)
N0 = 48;
fs0 = N0; % so f0=k gives exp(j*2*pi*k*n/N0)
n = (0:N0-1).';
K = 12;
U = zeros(N0, K);
for k = 1:K
    U(:, k) = sp_complex_tone(fs0, N0, k-1); % k-1 corresponds to DFT bin index
end

G = U' * U; % Gram matrix
maxOff = max(max(abs(G - diag(diag(G)))));
fprintf("orthogonality check: max off-diagonal |G| = %.3g\n", maxOff);
assert(maxOff < 1e-12, "Orthogonality check failed.");

figure;
imagesc(abs(G));
axis image;
colormap(gray);
cb = colorbar;
try
    cb.Color = [0 0 0];
catch
end
xlabel("k");
ylabel("k");
title("|U^H U| (orthogonality of complex exponentials)");

exportgraphics(gcf, fullfile(figDir, "ch03_orthogonality.png"), "Resolution", 200, "BackgroundColor", "white");

%% 6) angle and unwrap: phase "wrap" is a display convention
[x0, ~, nAx] = sp_complex_tone(fs, N, f0, "Phase", phi0); %#ok<ASGLU>
phi = angle(x0);
phiU = unwrap(phi);

Mshow = 400; % show first samples

figure;
tiledlayout(2, 1, "Padding", "compact", "TileSpacing", "compact");

nexttile;
plot(nAx(1:Mshow), phi(1:Mshow), "-", "LineWidth", 1.2);
grid on;
ylim([-pi, pi]);
yticks([-pi, -pi/2, 0, pi/2, pi]);
yticklabels({"-\pi", "-\pi/2", "0", "\pi/2", "\pi"});
xlabel("n");
ylabel("angle(x[n])");
title("angle: wraps to (-\pi,\pi]");

nexttile;
plot(nAx(1:Mshow), phiU(1:Mshow), "-", "LineWidth", 1.2);
hold on;
plot(nAx(1:Mshow), omega*nAx(1:Mshow) + phi0, "--", "LineWidth", 1.0);
grid on;
xlabel("n");
ylabel("unwrap(angle(x[n]))");
legend("unwrap", "\Omega n + \phi", "Location", "northwest");
title("unwrap: becomes (almost) a straight line whose slope is \Omega");

exportgraphics(gcf, fullfile(figDir, "ch03_phase_unwrap.png"), "Resolution", 200, "BackgroundColor", "white");
