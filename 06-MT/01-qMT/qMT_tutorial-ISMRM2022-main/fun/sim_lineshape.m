function [dataSim_lineshape, dataRaw_lineshape] = sim_lineshape()
  
lineshape = {'SuperLorentzian','Lorentzian','Gaussian'};

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
Opt.Method = 'Analytical equation';
Opt.ResetMz = false;

%% Varying lineshape

Opt.SNR = 1000;
for i=1:length(lineshape)
    Model.options.Lineshape = lineshape{i};
    Smodel = equation(Model, x, Opt);
    data.MTdata = addNoise(Smodel, Opt.SNR, 'mt');
    [SimCurveResults,Protocol,data] = getSimResults(Model, x, Opt, data);
    dataSim_lineshape(:,1,i) = SimCurveResults.Offsets;
    dataSim_lineshape(:,2,i) = SimCurveResults.curve(:,1);
    dataSim_lineshape(:,3,i) = SimCurveResults.curve(:,2);
    dataRaw_lineshape(:,1,i) = Protocol.Offsets;
    dataRaw_lineshape(:,2,i) = data.MTdata;
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