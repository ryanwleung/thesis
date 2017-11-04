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

for iFile = 1:nFile
    tempMatrix = xlsread(allNames{iFile});

    if size(tempMatrix, 2) ~= 4
        disp(allNames{iFile});
        error('Excel file does not have 4 columns of data');
    end
    
    MICP = table();
    MICP.SNW = tempMatrix(:, 3);
    MICP.PcGW = tempMatrix(:, 1);
    
    
    
    
%     tempTable = readtable(allNames{iFile});
%     tempDescriptionLogical = ~strncmp(tempTable.Properties.VariableNames, 'Var', 3);
%     tempTable.Properties.Description = tempTable.Properties.VariableNames{tempDescriptionLogical};
%     

%     for iColumn = 1:nColumn
%         while isempty( tempTable{1, iColumn} )
%             tempTable(1, iColumn) = [];
%         end
%     end
    
    
    
    
end














