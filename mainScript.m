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

test = 2;

if test == 1
    baseFlag = true;
    [ WaveParameterSensitivity , ~ , WaveBaseParameterSensitivity , dataBase ] = obj.RunSeismicAnalysisRoutine('ParameterSensitivity', baseFlag);
    
    obj.PlotParameterSensitivity(WaveParameterSensitivity, WaveBaseParameterSensitivity)
    obj.PlotPhaseSaturations()
    obj.PlotBackgroundProperties(dataBase)
    obj.PlotThicknessVsQuantity(WaveParameterSensitivity)
    obj.PlotPeakAmplitudeRatio(WaveParameterSensitivity, 'ParameterSensitivity')
    
elseif test == 2
    baseFlag = true;
    [ WaveOriginalResolution , data , WaveBaseOriginalResolution , ~ ] = obj.RunSeismicAnalysisRoutine('OriginalResolution', baseFlag);

    obj.PlotSeismogramOriginalResolution(WaveOriginalResolution)
    obj.PlotVelocityStructureOriginalResolution(WaveBaseOriginalResolution, WaveOriginalResolution)
    obj.PlotPeakAmplitudeRatio(WaveOriginalResolution, 'OriginalResolution')
    
end


