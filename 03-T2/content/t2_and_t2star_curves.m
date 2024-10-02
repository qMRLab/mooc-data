%% T2 and T2* decay curves
% Script to display T2 and T2* relaxometry curves for different tissues

% Simulation parameters
params.TE = linspace(0, 300, 100); % Echo times (in ms)

% Define T2 values for different tissues
params.T2_WM = 109.77; % T2 of white matter (in ms)
params.T2_GM = 96.07; % T2 of gray matter (in ms)

% Define T2* values for different tissues
params.T2star_WM = 67.63; % T2* of white matter (in ms)
params.T2star_GM = 48.48; % T2* of gray matter (in ms)

% Generate T2 and T2* decay signals
signal_WM_T2 = exp(-params.TE / params.T2_WM);
signal_GM_T2 = exp(-params.TE / params.T2_GM);
signal_WM_T2star = exp(-params.TE / params.T2star_WM);
signal_GM_T2star = exp(-params.TE / params.T2star_GM);

% Plot the T2 and T2* signals
figure;
hold on;
plot(params.TE, signal_WM_T2, '-b', 'DisplayName', 'T2 = 109.77 ms (white matter)');
plot(params.TE, signal_GM_T2, '-r', 'DisplayName', 'T2 = 96.07 ms (gray matter)');
plot(params.TE, signal_WM_T2star, '--b', 'DisplayName', 'T2* = 67.63 ms (white matter)');
plot(params.TE, signal_GM_T2star, '--r', 'DisplayName', 'T2* = 48.48 ms (gray matter)');
xlabel('Echo Time - TE (ms)');
ylabel('Transverse Magnetization (Mxy)');
legend();
title('T2 and T2* Decay Signals');

%% Export

TE = squeeze(params.TE)
save("t2_and_t2star_curvs.mat", "signal_WM_T2", "signal_GM_T2", "signal_WM_T2star", "signal_GM_T2star", "params", "EchoTimes")
