clc
close all
clear variables

load('180202 simulation results.mat')
figure
hold on
plot(seafloorDepthArray, minQuantityToFracture2PList, 'linewidth', 3)
plot(seafloorDepthArray, minQuantityToFracture3PList, 'linewidth', 3)
xlabel('Seafloor depth (m)')
ylabel('Minimum methane quantity to initiate a fracture (kg/m^3)')
legend('2P case')
legend('3P case')

figure
hold on
plot(seafloorDepthArray(1:end - 1), minQuantityToFracture2PList(2:end) - minQuantityToFracture2PList(1:end - 1), 'linewidth', 3)
plot(seafloorDepthArray(1:end - 1), minQuantityToFracture3PList(2:end) - minQuantityToFracture3PList(1:end - 1), 'linewidth', 3)
xlabel('Seafloor depth (m)')
ylabel('Difference in methane quantity to initiate a fracture (kg/m^3)')
legend('2P case')
legend('3P case')


% initialMethaneQuantity = 7;
% % initialMethaneQuantity = 122; % max methane quantity to still fracture
% 
% seafloorDepthArray = (500:50:2150)';
% % seafloorDepthArray = (2100:10:2200)';
% 
% %%% max seafloor depth for current parameters is 2150
% tic
% [minQuantityToFracture3PList, minQuantityToFracture2PList, errorList] = ...
%     DCTheoreticalFormation.RunMethaneQuantityFractureRoutine(seafloorDepthArray, initialMethaneQuantity);
% toc


% obj = DCTheoreticalFormation(2150, 0.4);
% [exportTable, transitionZoneProperties] = obj.RunSolubilitySaturationRoutine(122);
% exportTable = obj.RunRockAndRatioRoutine(exportTable);
% obj.GenerateResultPlots(exportTable, transitionZoneProperties);








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

