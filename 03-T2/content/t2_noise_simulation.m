%% T2 and T2* decay curves

% Define model
Model = mono_t2;

% Define signal parameters for different tissues
x = struct;
x.M0 = 1000;
x.T2 = 109; % (in ms)

% Define the signal-to-noise ratio 
Opt1.SNR = 10;
Opt2.SNR = 50;
Opt3.SNR = 90;
Opt4.SNR = 130;

% Run the simulation for T2 and T2* decay curves
[FitResult_SNR10, data_SNR10] = Model.Sim_Single_Voxel_Curve(x, Opt1);
[FitResult_SNR50, data_SNR50] = Model.Sim_Single_Voxel_Curve(x, Opt2);
[FitResult_SNR90, data_SNR90] = Model.Sim_Single_Voxel_Curve(x, Opt3);
[FitResult_SNR130, data_SNR130] = Model.Sim_Single_Voxel_Curve(x, Opt4);

% T2 constants
T2_SNR10 = FitResult_SNR10.T2;
T2_SNR50 = FitResult_SNR50.T2;
T2_SNR90 = FitResult_SNR90.T2;
T2_SNR130 = FitResult_SNR130.T2;

% T2 decay curves
signal_SNR10 = FitResult_SNR10.M0/1000 * exp(-params.TE / FitResult_SNR10.T2);
signal_SNR50 = FitResult_SNR50.M0/1000 * exp(-params.TE / FitResult_SNR50.T2);
signal_SNR90 = FitResult_SNR90.M0/1000 * exp(-params.TE / FitResult_SNR90.T2);
signal_SNR130 = FitResult_SNR130.M0/1000 * exp(-params.TE / FitResult_SNR130.T2);

% Noisy data points
EchoTimes  = [12.8000; 25.6000; 38.4000; 51.2000; 64.0000; 76.8000; 89.6000; 102.4000; 115.2000; 128.0000; 140.8000; 153.6000; 166.4000; 179.2000; 192.0000; 204.8000; 217.6000; 230.4000; 243.2000; 256.0000; 268.8000; 281.6000; 294.4000; 307.2000; 320.0000; 332.8000; 345.6000; 358.4000; 371.2000; 384.0000];
SEdata_SNR10 = data_SNR10.SEdata/1000;
SEdata_SNR50 = data_SNR50.SEdata/1000;
SEdata_SNR90 = data_SNR90.SEdata/1000;
SEdata_SNR130 = data_SNR130.SEdata/1000;

%% Export

disp(EchoTimes)
disp(params.TE)

save("t2_noise_simulation.mat", "signal_SNR10", "signal_SNR50", "signal_SNR90", "signal_SNR130", "T2_SNR10", "T2_SNR50", "T2_SNR90", "T2_SNR130", "SEdata_SNR10", "SEdata_SNR50", "SEdata_SNR90", "SEdata_SNR130", "EchoTimes", "params")
