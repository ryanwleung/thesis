clc
close all
clear variables

ch4Quantity = 40;

% formationInstance = DCHydrateRidge();
formationInstance = DCBlakeRidge();

[ exportTable , transitionZoneProperties ] = formationInstance.RunSolubilitySaturationRoutine( ch4Quantity );

formationInstance.GenerateResultPlots( exportTable , transitionZoneProperties );
formationInstance.PlotMICP();
formationInstance.PlotCumPSD();
formationInstance.PlotPSD( 'linear' );
formationInstance.PlotPSD( 'log' );



