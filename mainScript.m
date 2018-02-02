clc
close all
clear variables



n = 2;
seafloorDepthArray = linspace(1500, 1600, n);

[minQuantityToFracture3PList, minQuantityToFracture2PList, errorList] = DCTheoreticalFormation.RunMethaneQuantityFractureRoutine(seafloorDepthArray);

% ch4Quantity = 115; % near max for KB
% ch4Quantity = 55; % near max for HR
% ch4Quantity = 120; % near max for BR
% ch4Quantity = 40;

% obj = DCHydrateRidge();
% obj = DCBlakeRidge();
% obj = DCKumanoBasin();

% [exportTable, transitionZoneProperties] = obj.RunSolubilitySaturationRoutine(ch4Quantity);
% exportTable = obj.RunRockAndRatioRoutine(exportTable);

% obj.GenerateResultPlots(exportTable, transitionZoneProperties);


%%% Plot results
% obj.PlotMICP();
% obj.PlotCumPSD();
% obj.PlotPSD('linear');
% obj.PlotPSD('log');
% obj.CalcPoreVolumeDistribution();


%%% Utility calls
% % Only run this static method when you need to update the way
% % the MICP is saved in MATLAB
% DCKumanoBasin.ExtractMICP();






% %%% Seismic analysis for Blake Ridge
% obj = DCSeismicAnalysisBR();
% 
% test = 1;
% 
% switch test
%     case 1
%         [ WaveParameterSensitivity , ~ , WaveBaseParameterSensitivity , dataBase , ~ ] = obj.RunSeismicAnalysisRoutine('ParameterSensitivity');
% 
%         obj.PlotParameterSensitivity(WaveParameterSensitivity, WaveBaseParameterSensitivity)
%         obj.PlotPhaseSaturations()
%         obj.PlotBackgroundProperties(dataBase)
%         obj.PlotThicknessVsQuantity(WaveParameterSensitivity)
% %         obj.PlotPeakAmplitudeRatio(WaveParameterSensitivity)
%         
%     case 2
%         obj.Dickens = obj.LoadDickensBlakeRidge();
%         [ WaveOriginalResolution , data , WaveBaseOriginalResolution , ~ , WaveDickens ] = obj.RunSeismicAnalysisRoutine('OriginalResolution');
%         
%         obj.PlotSeismogramOriginalResolution(WaveOriginalResolution)
%         obj.PlotVelocityStructureOriginalResolution(WaveBaseOriginalResolution, WaveOriginalResolution)
% %         obj.PlotPeakAmplitudeRatio(WaveOriginalResolution)
%         obj.PlotDickensSeismogram(WaveOriginalResolution, WaveDickens)
%         
% end

