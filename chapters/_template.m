%% Chapter Template
% このファイルは chapters/ 配下の章スクリプトの雛形。
%
% - 章スクリプトは「Run All」で図まで再現できる形にする
% - 必要なときだけ rng(0) を使う

clear; close all; clc;

% ルート（bookRoot）を推定して src/ を path に追加
thisFile = mfilename("fullpath");
bookRoot = fileparts(fileparts(thisFile));
addpath(fullfile(bookRoot, "src"));

% 出力先（図）
figDir = fullfile(bookRoot, "figs");
if ~exist(figDir, "dir")
    mkdir(figDir);
end

%% Parameters (example)
fs = 8000;   % [Hz]
N  = 1024;   % number of samples

%% Experiment (TODO)
% TODO: write experiment here

%% Checks (TODO)
% TODO: add assert-based checks here
