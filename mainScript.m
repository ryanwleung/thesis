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
obj.RunSeismicAnalysisRoutine();













