clc
close all
clear variables

% load('180202 simulation results.mat')
% load('180211 simulation results.mat')
% load('180219 simulation results.mat')
% load('180307 simulation results.mat')
% DCTheoreticalFormation.PlotCombined(seafloorDepthArray, ...
%                                     minQuantityToFracture2PList, minQuantityToFracture3PList, ...
%                                     sgFracture2PList, sgFracture3PList)

% DCTheoreticalFormation.PlotMethaneQuantities(seafloorDepthArray, minQuantityToFracture2PList, minQuantityToFracture3PList);
% DCTheoreticalFormation.PlotGasSaturations(seafloorDepthArray, sgFracture2PList, sgFracture3PList);
% DCTheoreticalFormation.PlotDepths(seafloorDepthArray, depthStructList);

% PlotAllMICP()



%%% Daigle scenario generation
obj = DCKumanoBasin();
obj.runDaigleCases = true;
obj.scenario = 4;

switch obj.scenario
    case 1
        ch4Quantity = 100;
    case 2
        ch4Quantity = 27;
    case 3
        ch4Quantity = 53;
    case 4
        ch4Quantity = 38;
end

[exportTable, transitionZoneProperties] = obj.RunSolubilitySaturationRoutineExperimental(ch4Quantity);
obj.PlotDaigleScenario(exportTable);



%%% max seafloor depth for current parameters is 2100
% seafloorDepthArray = (500:50:2100)';
% seafloorDepthArray = (500:50:800)';
% initialMethaneQuantity = 7;


% seafloorDepthArray = 600;
% initialMethaneQuantity = 7.6; % max methane quantity to still fracture

% tic
% [minQuantityToFracture3PList, minQuantityToFracture2PList, sgFracture3PList, sgFracture2PList, errorList, depthStructList] = ...
%     DCTheoreticalFormation.RunMethaneQuantityFractureRoutine(seafloorDepthArray, initialMethaneQuantity);
% toc

% tempDepth = [500; 550; 600];
% tempMethaneQuantity = [7.4];
% 
% 
% obj = DCTheoreticalFormation(550, 0.4);
% [exportTable, transitionZoneProperties] = obj.RunSolubilitySaturationRoutine(7.4);

% obj = DCTheoreticalFormation(1800, 0.4);
% [exportTable, transitionZoneProperties] = obj.RunSolubilitySaturationRoutine(60);
% 
% exportTable = obj.RunRockAndRatioRoutine(exportTable);
% obj.GenerateResultPlots(exportTable, transitionZoneProperties);








% ch4Quantity = 115; % near max for KB
% ch4Quantity = 55; % near max for HR
% ch4Quantity = 120; % near max for BR

% ch4Quantity = 40;
% % 
% % obj = DCBlakeRidge();
% % obj = DCHydrateRidge();
% obj = DCKumanoBasin();
% 
% [exportTable, transitionZoneProperties] = obj.RunSolubilitySaturationRoutine(ch4Quantity);
% exportTable = obj.RunRockAndRatioRoutine(exportTable);
% 
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
% %         obj.PlotPhaseSaturations()
%         obj.PlotBackgroundProperties(dataBase)
%         obj.PlotThicknessVsQuantity(WaveParameterSensitivity)
% %         obj.PlotPeakAmplitudeRatio(WaveParameterSensitivity)
%         obj.PlotTWTTLength(WaveParameterSensitivity)
%         
%     case 2
%         obj.Dickens = obj.LoadDickensBlakeRidge();
%         [ WaveOriginalResolution , data , WaveBaseOriginalResolution , ~ , WaveDickens ] = obj.RunSeismicAnalysisRoutine('OriginalResolution');
%         
%         obj.PlotSeismogramOriginalResolution(WaveOriginalResolution)
%         obj.PlotVelocityStructureOriginalResolution(WaveBaseOriginalResolution, WaveOriginalResolution)
% %         obj.PlotPeakAmplitudeRatio(WaveOriginalResolution)
%         obj.PlotTWTTLength(WaveOriginalResolution)
%         obj.PlotDickensSeismogram(WaveOriginalResolution, WaveDickens)
%         
% end

function PlotAllMICP()
    figure
    
    

    obj = DCBlakeRidge();
    semilogy(1 - obj.MICP1.S_nw, obj.MICP1.Pc_gw, 'Linewidth', 4)
    hold on
    semilogy(1 - obj.MICP2.S_nw, obj.MICP2.Pc_gw, 'Linewidth', 4)
    
    obj = DCHydrateRidge();
    % Copied from DCHydrateRidge method
    semilogy(1 - obj.MICP1.S_nw, obj.MICP1.Pc_gw, 'Linewidth', 4)
    
    
    
    obj = DCKumanoBasin();
    semilogy(1 - obj.MICP{1}.SNW, obj.MICP{1}.PcGW, 'Linewidth', 4)
    
%     xlabel('1 - S_n_w, or S_w')
    xlabel('S_w')
    ylabel('P_c_g_w (MPa)')
%     title('Primary Drainage Capillary Pressure')
    legend('Blake Ridge 1', 'Blake Ridge 2', 'Hydrate Ridge', 'Kumano Basin')
end