% Define the model 
Model = mwf;

% Load data into environment, and rotate mask to be aligned with IR data
data = struct;
load('../data/mwf/MET2data.mat');
load('../data/mwf/Mask.mat');
data.MET2data = double(MET2data);
data.Mask = double(Mask);

% Define fitting parameters
EchoTime  = [12.8000; 25.6000; 38.4000; 51.2000; 64.0000; 76.8000; 89.6000; 102.4000; 115.2000; 128.0000; 140.8000; 153.6000; 166.4000; 179.2000; 192.0000; 204.8000; 217.6000; 230.4000; 243.2000; 256.0000; 268.8000; 281.6000; 294.4000; 307.2000; 320.0000; 332.8000; 345.6000; 358.4000; 371.2000; 384.0000];
Model.Prot.SEdata.Mat = [ EchoTime ];

% MET2w MRI data at different TE values
ME_TE_1 = imrotate(squeeze(data.MET2data(:,:,:,1).*data.Mask),-90);
ME_TE_2 = imrotate(squeeze(data.MET2data(:,:,:,10).*data.Mask),-90);
ME_TE_3 = imrotate(squeeze(data.MET2data(:,:,:,20).*data.Mask),-90);
ME_TE_4 = imrotate(squeeze(data.MET2data(:,:,:,30).*data.Mask),-90);

% Fit the data
FitResults_mwf = FitData(data,Model,0);

MWF = imrotate(squeeze(FitResults_mwf.MWF.*data.Mask), -90);
T2MW = imrotate(squeeze(FitResults_mwf.T2MW.*data.Mask), -90);
T2IEW = imrotate(squeeze(FitResults_mwf.T2IEW.*data.Mask), -90);

%% Export

save("multiexpo_T2_image.mat", "ME_TE_1", "ME_TE_2", "ME_TE_3", "ME_TE_4", "EchoTime", "FitResults_mwf", "MWF", "T2MW", "T2IEW")
