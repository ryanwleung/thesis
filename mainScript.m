clc
close all
clear variables

% ch4Quantity = 40;
% 
% % formationInstance = DCHydrateRidge();
% formationInstance = DCBlakeRidge();
% 
% [ exportTable , transitionZoneProperties ] = formationInstance.RunSolubilitySaturationRoutine( ch4Quantity );
% 
% formationInstance.GenerateResultPlots( exportTable , transitionZoneProperties );
% formationInstance.PlotMICP();
% formationInstance.PlotCumPSD();
% formationInstance.PlotPSD( 'linear' );
% formationInstance.PlotPSD( 'log' );


%%% Only run this static method when you need to update the way
%%% the MICP is saved in MATLAB
% DCKumanoBasin.ExtractMICP();

% formationInstance = DCKumanoBasin();
% formationInstance.PlotMICP();
% formationInstance.PlotCumPSD();
% % formationInstance.PlotPSD( 'linear' );
% formationInstance.PlotPSD( 'log' );


obj = DCSeismicAnalysisBR();



baseFlag = true;
[ WaveParameterSensitivity , ~ , WaveBase , dataBase ] = obj.RunSeismicAnalysisRoutine('ParameterSensitivity', baseFlag);

obj.PlotParameterSensitivity(WaveParameterSensitivity, WaveBase)
% obj.PlotPhaseSaturations(WaveParameterSensitivity)
% 
% obj.PlotBackgroundProperties(dataBase)
% 
% obj.PlotThicknessVsQuantity(WaveParameterSensitivity)

% baseFlag = false;
% [ WaveOriginalResolution , data , ~ , ~ ] = obj.RunSeismicAnalysisRoutine('OriginalResolution', baseFlag);
% 
% obj.PlotSeismogramOriginalResolution(WaveOriginalResolution)



