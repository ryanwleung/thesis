clc
close all
clear variables

allFiles = dir;
allNames = { allFiles(~[allFiles.isdir]).name }';

allNamesLogical = contains(allNames, '.XLSX');
allNames = allNames(allNamesLogical);

nFile = numel(allNames);
if nFile ~= 19
    error('Number of Excel files does not match 19')
end

MICPCellArray = cell(nFile, 1);



sigmaHgAir = 485; % dynes/cm
sigmaGasWater = 72; % dynes/cm
sigmaHydrateWater = 27; % dynes/cm

thetaHgAir = 140; % deg
thetaGasWater = 180; % deg
thetaHydrateWater = 180; % deg

conversionPsiToMpa = 1/145.03774; % multiply by this

for iFile = 1:nFile
    tempMatrix = xlsread(allNames{iFile});

    if size(tempMatrix, 2) ~= 4
        disp(allNames{iFile});
        error('Excel file does not have 4 columns of data');
    end
    
    MICP = table();
    MICP.SNW = tempMatrix(:, 2);
    
    MICP.PcGW = tempMatrix(:, 1) .* (sigmaGasWater * cosd(thetaGasWater)) ...
                                  / (sigmaHgAir * cosd(thetaHgAir)); % psi
    MICP.PcHW = tempMatrix(:, 1) .* (sigmaHydrateWater * cosd(thetaHydrateWater)) ...
                                  / (sigmaHgAir * cosd(thetaHgAir)); % psi
    
    MICP.PcGW = MICP.PcGW .* conversionPsiToMpa;
    MICP.PcHW = MICP.PcHW .* conversionPsiToMpa;
    
    % Convert pore throat radius in microns to pore throat diameter in meters
    MICP.PoreThroatDiameter = tempMatrix(:, 4) .* 2 ./ 1e6;
    
    MICPCellArray{iFile} = MICP;
end














