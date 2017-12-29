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

test = 1;

switch test
    case 1
        [ WaveParameterSensitivity , ~ , WaveBaseParameterSensitivity , dataBase , ~ ] = obj.RunSeismicAnalysisRoutine('ParameterSensitivity');

        obj.PlotParameterSensitivity(WaveParameterSensitivity, WaveBaseParameterSensitivity)
        obj.PlotPhaseSaturations()
        obj.PlotBackgroundProperties(dataBase)
        obj.PlotThicknessVsQuantity(WaveParameterSensitivity)
        obj.PlotPeakAmplitudeRatio(WaveParameterSensitivity)
    
    case 2
        obj.Dickens = obj.LoadDickensBlakeRidge();
        [ WaveOriginalResolution , data , WaveBaseOriginalResolution , ~ , WaveDickens ] = obj.RunSeismicAnalysisRoutine('OriginalResolution');
        
        obj.PlotSeismogramOriginalResolution(WaveOriginalResolution)
        obj.PlotVelocityStructureOriginalResolution(WaveBaseOriginalResolution, WaveOriginalResolution)
        obj.PlotPeakAmplitudeRatio(WaveOriginalResolution)
        obj.PlotDickensSeismogram(WaveOriginalResolution, WaveDickens)
        
end


