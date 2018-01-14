clc
close all
clear variables

ch4Quantity = 40;

obj = DCHydrateRidge();
% obj = DCBlakeRidge();
% obj = DCKumanoBasin();


[exportTable, transitionZoneProperties] = obj.RunSolubilitySaturationRoutine(ch4Quantity);


%%% Plot results
% obj.PlotMICP();
% obj.PlotCumPSD();
% obj.PlotPSD('linear');
% obj.PlotPSD('log');
% obj.CalcPoreVolumeDistribution();

obj.GenerateResultPlots(exportTable, transitionZoneProperties);

%%% Utility calls
% % Only run this static method when you need to update the way
% % the MICP is saved in MATLAB
% DCKumanoBasin.ExtractMICP();


%%% Seismic analysis for Blake Ridge
% obj = DCSeismicAnalysisBR();
% 
% test = 2;
% 
% switch test
%     case 1
%         [ WaveParameterSensitivity , ~ , WaveBaseParameterSensitivity , dataBase , ~ ] = obj.RunSeismicAnalysisRoutine('ParameterSensitivity');
% 
%         obj.PlotParameterSensitivity(WaveParameterSensitivity, WaveBaseParameterSensitivity)
%         obj.PlotPhaseSaturations()
%         obj.PlotBackgroundProperties(dataBase)
%         obj.PlotThicknessVsQuantity(WaveParameterSensitivity)
%         obj.PlotPeakAmplitudeRatio(WaveParameterSensitivity)
%         
%     case 2
%         obj.Dickens = obj.LoadDickensBlakeRidge();
%         [ WaveOriginalResolution , data , WaveBaseOriginalResolution , ~ , WaveDickens ] = obj.RunSeismicAnalysisRoutine('OriginalResolution');
%         
%         obj.PlotSeismogramOriginalResolution(WaveOriginalResolution)
%         obj.PlotVelocityStructureOriginalResolution(WaveBaseOriginalResolution, WaveOriginalResolution)
%         obj.PlotPeakAmplitudeRatio(WaveOriginalResolution)
%         obj.PlotDickensSeismogram(WaveOriginalResolution, WaveDickens)
%         
% end

% a = linspace(0.0022, 0, 10)';
% b = zeros(numel(a), 1);
% c = zeros(numel(a), 1);
% radiusG = 9.307991980734942e-08;
% radiusH = 9.558501814783801e-08;
% for i = 1:numel(a)
%     [ b(i) , c(i) ] = obj.InterpCumPSD( a(i) , radiusG , radiusH );
% end