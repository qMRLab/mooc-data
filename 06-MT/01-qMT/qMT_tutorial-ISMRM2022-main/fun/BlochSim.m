function [datablochSim_Pulses, dataSim_Pulses, dataRaw_Pulses] = BlochSim(numPulses)
  
% Define qMT Model
Model = qmt_spgr;

% Input parameters
x = struct;
x.F = 0.16;
x.kr = 30;
x.R1f = 1;
x.R1r = 1;
x.T2f = 0.03;
x.T2r = 1.3e-05;

% Set simulation options
Opt.SNR = 1000;
Opt.Method = 'Block equation';
Opt.ResetMz = false;

%% Varying number of pulses (Bloch Simulations)
datablochSim_Pulses = zeros(5,3,length(numPulses));
for i=1:length(numPulses)
    Model.options.MT_Pulse_NofMTpulses = numPulses(i);
    %Bloch simulation data
    Smodel = equation(Model, x, Opt);
    blochSim_angle1 = Smodel(1:2:end,1);
    blochSim_angle2 = Smodel(2:2:end,1);
    data.MTdata = addNoise(Smodel, Opt.SNR, 'mt');
    %Analytical solution data
    [SimCurveResults,Protocol,data] = getSimResults(Model, x, Opt, data);
    
    dataSim_Pulses(:,1,i) = SimCurveResults.Offsets;
    dataSim_Pulses(:,2,i) = SimCurveResults.curve(:,1);
    dataSim_Pulses(:,3,i) = SimCurveResults.curve(:,2);
    datablochSim_Pulses(:,1,i) = unique(Protocol.Offsets);
    datablochSim_Pulses(:,2,i) = blochSim_angle1;
    datablochSim_Pulses(:,3,i) = blochSim_angle2;
    dataRaw_Pulses(:,1,i) = Protocol.Offsets;
    dataRaw_Pulses(:,2,i) = data.MTdata;
end
end

%% Function to get the results from the simulations varying different parameters
function [SimCurveResults,Protocol,data] = getSimResults(Model, x, Opt, data)


    Protocol = GetProt(Model);
    FitOpt   = GetFitOpt(Model,data);
    % normalize data
            NoMT = Protocol.Angles<1;
            if ~any(NoMT)
                warning('No MToff (i.e. no volumes acquired with Angles=0) --> Fitting assumes that MTData are already normalized.');
            else
                data.MTdata = data.MTdata/median(data.MTdata(NoMT));
                data.MTdata = data.MTdata(~NoMT);
                Protocol.Angles  = Protocol.Angles(~NoMT)*0.5;
                Protocol.Offsets = Protocol.Offsets(~NoMT);
            end
            SimCurveResults = SPGR_SimCurve(x, Protocol, FitOpt );
end
