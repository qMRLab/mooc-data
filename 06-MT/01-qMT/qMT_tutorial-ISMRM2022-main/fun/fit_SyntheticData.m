function [b0FittedResults,b0NormFittedResults,b1FittedResults,b1NormFittedResults,t1FittedResults,t1NormFittedResults] = fit_SyntheticData(B0map,B1map,T1map)

R1map = 1./T1map;

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

%% Fitted parameters
Opt.SNR = 1000;
Smodel = equation(Model, x, Opt);
data.MTdata = addNoise(Smodel, Opt.SNR, 'mt');
[SimCurveResults,Protocol,data] = getSimResults(Model, x, Opt, data);

mtdata = reshape(data.MTdata,[1,1,1,10]);

B1map_default = 1;
B0map_default = 0;
R1map_default = 1/convert_T1f_T1meas([x.F, x.kr*x.F, 1/x.R1f, 1/x.R1r],'t1f_2_t1meas');

%Vary B1map
for i=1:length(B1map)
    data.B1map(1,1,1) = B1map(i);
    data.B0map(1,1,1) = B0map_default;
    data.R1map(1,1,1) = R1map_default;
    data.MTdata = mtdata;
    
    FitResults = FitData(data,Model,0);
    
    fn = fieldnames(FitResults);
    for m=1:length(fn)
        B1FitResults{i,m} = FitResults.(fn{m});
    end
end

%Vary B0map
for i=1:length(B0map)
    data.B1map(1,1,1) = B1map_default;
    data.B0map(1,1,1) = B0map(i);
    data.R1map(1,1,1) = R1map_default;
    data.MTdata = mtdata;
    
    FitResults = FitData(data,Model,0);
    
    fn = fieldnames(FitResults);
    for m=1:length(fn)
        B0FitResults{i,m} = FitResults.(fn{m});
    end
end

%Vary R1map
for i=1:length(R1map)
    data.B1map(1,1,1) = B1map_default;
    data.B0map(1,1,1) = B0map_default;
    data.R1map(1,1,1) = R1map(i);
    data.MTdata = mtdata;
    
    FitResults = FitData(data,Model,0);
    
    fn = fieldnames(FitResults);
    for m=1:length(fn)
        T1FitResults{i,m} = FitResults.(fn{m});
    end
end

%Get fitted parameters
b1FittedResults = horzcat(B1map',[B1FitResults{:,1}]',[B1FitResults{:,2}]',[B1FitResults{:,3}]',[B1FitResults{:,4}]',[B1FitResults{:,5}]',[B1FitResults{:,6}]',[B1FitResults{:,7}]',[B1FitResults{:,8}]');
b0FittedResults = horzcat(B0map',[B0FitResults{:,1}]',[B0FitResults{:,2}]',[B0FitResults{:,3}]',[B0FitResults{:,4}]',[B0FitResults{:,5}]',[B0FitResults{:,6}]',[B0FitResults{:,7}]',[B0FitResults{:,8}]');
t1FittedResults = horzcat(1./R1map',[T1FitResults{:,1}]',[T1FitResults{:,2}]',[T1FitResults{:,3}]',[T1FitResults{:,4}]',[T1FitResults{:,5}]',[T1FitResults{:,6}]',[T1FitResults{:,7}]',[T1FitResults{:,8}]');

b1NormFittedResults(:,1) = B1map';
b1NormFittedResults(:,2) = (b1FittedResults(:,2) - x.F)*100 / (x.F);
b1NormFittedResults(:,3) = (b1FittedResults(:,8) - x.kr*x.F)*100 / (x.kr*x.F);
b1NormFittedResults(:,4) = (b1FittedResults(:,6) - x.T2f)*100 / (x.T2f);
b1NormFittedResults(:,5) = (b1FittedResults(:,7) - x.T2r)*100 / (x.T2r);

b0NormFittedResults(:,1) = B0map';
b0NormFittedResults(:,2) = (b0FittedResults(:,2) - x.F)*100 / (x.F);
b0NormFittedResults(:,3) = (b0FittedResults(:,8) - x.kr*x.F)*100 / (x.kr*x.F);
b0NormFittedResults(:,4) = (b0FittedResults(:,6) - x.T2f)*100 / (x.T2f);
b0NormFittedResults(:,5) = (b0FittedResults(:,7) - x.T2r)*100 / (x.T2r);

t1NormFittedResults(:,1) = 1./R1map';
t1NormFittedResults(:,2) = (t1FittedResults(:,2) - x.F)*100 / (x.F);
t1NormFittedResults(:,3) = (t1FittedResults(:,8) - x.kr*x.F)*100 / (x.kr*x.F);
t1NormFittedResults(:,4) = (t1FittedResults(:,6) - x.T2f)*100 / (x.T2f);
t1NormFittedResults(:,5) = (t1FittedResults(:,7) - x.T2r)*100 / (x.T2r);

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