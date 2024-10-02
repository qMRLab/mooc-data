Model = mwf

% Define initial MWF and T2 times of myelin water and intra- and extracellular water
x = struct;
x.MWF = 50;
x.T2MW = 20;
x.T2IEW = 120;


% Define echo times
params.TE = linspace(0, 300, 100);

% Set simulation options
Opt.SNR = 120;
Opt.T2Spectrumvariance_Myelin = 5;
Opt.T2Spectrumvariance_IEIntraExtracellularWater = 20;

% Run simulation
figure('Name','Single Voxel Curve Simulation');
FitResult = Model.Sim_Single_Voxel_Curve(x,Opt);

% T2 relaxation curves for myelin water and intra/extracellular water
% (using a mono-exponential curve)
signal_mono_MW = exp(-params.TE / FitResult.T2MW);
signal_mono_IEW = exp(-params.TE / FitResult.T2IEW);

% T2 relaxation curve for multi-expo model
signal_multi_MWF = (FitResult.MWF/100)*signal_mono_MW + (1 - FitResult.MWF/100)*signal_mono_IEW;

%% Export

TE = squeeze(params.TE)
save("multiexpo_T2_curves.mat", "signal_mono_MW", "signal_mono_IEW", "TE", "FitResult", "params", "x", "Opt")
