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
<<<<<<< HEAD
% formationInstance.
PlotCumPSD();
% % formationInstance.PlotPSD( 'linear' );
% formationInstance.PlotPSD( 'log' );
<<<<<<< HEAD


seismicInstance = DCSeismicAnalysisBR();
seismicInstance.RunSeismicAnalysisRoutine();








=======
>>>>>>> 75fb1eba223d99279a20e43402f14bd6462377f2
=======
% formationInstance.PlotCumPSD();
% % formationInstance.PlotPSD( 'linear' );
% formationInstance.PlotPSD( 'log' );
>>>>>>> 75fb1eba223d99279a20e43402f14bd6462377f2


obj = DCSeismicAnalysisBR();
obj.RunSeismicAnalysisRoutine();













