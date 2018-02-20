classdef DCTheoreticalFormation < BCFormation
    properties
        MICP % Cell array of tables
        colorOrder % array containing order index of loaded MICP tables
        depthOrder % array containing depths in the order of loaded MICP tables
        MICPInterp % Cell array of tables of MICP that will actually be used
        
        poissonsRatio
        solAxis
    end
    properties (Constant)
%         satAxis = [0 1 0 750];
%         pcgwAxis = [0 2.5 0 750];
%         ratioAxis = [0 1.2 0 750];
        
        satAxis = [0 1 365 400];
        pcgwAxis = [0 2.5 365 400];
        ratioAxis = [0 1.2 365 400];
        solAxisTemp = [0.153 0.166 365 400];
        
    end
    methods
        %%% Constructor
        function [ obj ] = DCTheoreticalFormation(seafloorDepth_, poissonsRatio_)
            obj@BCFormation()
            
            obj.poissonsRatio = poissonsRatio_;
            
            
            obj.seafloorDepth = seafloorDepth_; % m
            obj.minDepth = 1;                   % mbsf
            obj.maxDepth = 750;                 % mbsf
            obj.depthIncrement = 0.5;           % m
            obj.depthArray = (obj.minDepth : obj.depthIncrement : obj.maxDepth)';
            
            obj.temperatureGradient = 40;       % C deg/km (C0002 - Daigle and Dugan 2014)            
            obj.seafloorTemperature = 3;    % C deg
            obj.salinityWtPercent = 3.5;        % weight percent (wt%) of NaCl in seawater
            
            obj.phi0 = 0.66;
            obj.phiInf = 0.1333;
            obj.B = 1333; % m
            
            
            
            
            
            % Loads MICP from a .mat file into the object property MICP
            obj.MICP = DCTheoreticalFormation.LoadMICP();
            
            
            
            obj.MICPInterp = DCTheoreticalFormation.SelectMICPForInterp(obj.MICP);
%             
%             % Testing plots only on selected MICP data for interpolation
%             obj.MICP = obj.MICPInterp;
        end
        
        %%% Petrophysical calculations
        function [ sg3PAtFracture , sg2PAtFracture ] = CheckFractureStatus( ~ , exportTable , transitionZoneProperties )
            sg3PAtFracture = [];
            sg2PAtFracture = [];
            
            iTopGas3P = find(exportTable.GasSat3P > 0, 1);
            iTopGas2P = transitionZoneProperties.Bulk3PSolEQLIndex;
            
            fractured3PIndex = find(exportTable.ratio3P(iTopGas3P:end) >= 1, 1);
            fractured2PIndex = find(exportTable.ratio2P(iTopGas2P:end) >= 1, 1);
            
            if ~isempty(fractured3PIndex)
                actual3PIndex = iTopGas3P + fractured3PIndex - 1;
                sg3PAtFracture = exportTable.GasSat3P(actual3PIndex);
            end
            if ~isempty(fractured2PIndex)
                actual2PIndex = iTopGas2P + fractured2PIndex - 1;
                sg2PAtFracture = exportTable.GasSat2P(actual2PIndex);
            end
        end
        function [ slope ] = CalcSlopeOfCumPSD( obj, stringType )
            % slope is in volume fraction per radius meter
            
            n = numel(obj.MICP);
            slope = cell(n, 1);
            
            for i = 1:n
                diameterInMeters = obj.MICP{i}.PoreThroatDiameter;
                radiusInMeters = diameterInMeters ./ 2;
                sNw = obj.MICP{i}.SNW;
                
                switch stringType
                    case 'linear'
                        slope{i} =  -( sNw(2:end) - sNw(1:end - 1) ) ...
                                  ./ ( radiusInMeters(2:end) - radiusInMeters(1:end - 1) );
                    case 'log'
                        slope{i} =  -( sNw(2:end) - sNw(1:end - 1) ) ...
                                  ./ ( log10(radiusInMeters(2:end)) - log10(radiusInMeters(1:end - 1)) );
                end
            end
        end
        
        
        %%% Plotting subclass functions
        function PlotSolLooped( obj , solFigure , exportTable , transitionZoneProperties , doPlotBulkAndMinSol )
            
            expandAxisLimitFactor = 0.1;
            iTop = transitionZoneProperties.Top3PIndex;
            iBottom = transitionZoneProperties.Bottom3PIndex;
            
            minSolubility = exportTable.OverallSol(iTop) * (1 - expandAxisLimitFactor);
            maxSolubility = exportTable.OverallSol(iBottom) * (1 + expandAxisLimitFactor);
            
            minDepth = exportTable.Depth(iTop) * (1 - expandAxisLimitFactor);
            maxDepth = exportTable.Depth(iBottom) * (1 + expandAxisLimitFactor);
            
            
            
            obj.solAxis = [minSolubility, maxSolubility, minDepth, maxDepth];
            obj.PlotSol( solFigure , exportTable , doPlotBulkAndMinSol );
        end
        function PlotMICP( obj )
            figure
            n = numel(obj.MICP);            
            for i = 1:n
                PcGW = obj.MICP{i}.PcGW;
                SNW = obj.MICP{i}.SNW;
                
                hold on
                plot(SNW, PcGW, 'Linewidth', 1)
                
                textTable = obj.MICP{i};
                textTable(textTable.SNW == 0, :) = [];
                textIndex = round(i / n * height(textTable));
                text(textTable.SNW(textIndex), textTable.PcGW(textIndex), num2str(obj.MICP{i}.Properties.UserData))
            end
            xlabel('1 - S_n_w, or S_w')
            ylabel('Pc in MPa')
            title('Primary Drainage Capillary Pressure Curve')
            set(gca, 'XDir', 'reverse')
            set(gca, 'Yscale', 'log')
            
            
            [averagePcgwPlot, snwPlot] = obj.GetInterpMICPForPlotting();
            scatter(snwPlot, averagePcgwPlot, 75)
            
        end
        function PlotCumPSD( obj )
            n = numel(obj.MICP);
            
            figure
            hold on
            for i = 1:n
                diameterInMeters = obj.MICP{i}.PoreThroatDiameter;
                radiusInMeters = diameterInMeters ./ 2;
                SNW = obj.MICP{i}.SNW;
                
                plot(radiusInMeters, SNW, 'Linewidth', 1)
                
                textTable = obj.MICP{i};
                textTable.PoreThroatDiameter = radiusInMeters; % hard fix for swapping to radius from diameter
                textTable(textTable.SNW == 0, :) = [];
                textIndex = round(i / n * height(textTable));
                text(textTable.PoreThroatDiameter(textIndex), textTable.SNW(textIndex), num2str(obj.MICP{i}.Properties.UserData))
            end
            xlabel('Pore radius in meters')
            ylabel('S_n_w')
            title('Cumulative Pore Size Distribution')
            axis([1e-9, 1e-4, 0, 1])
            set(gca, 'XDir', 'reverse')
            set(gca, 'Xscale', 'log')
        end
        function PlotPSD( obj , stringType )
            % slope is a cell array containing the double array of each
            % slope from each MICP data set
            slope = obj.CalcSlopeOfCumPSD(stringType);
            n = numel(obj.MICP);
            
            figure
            hold on
            for i = 1:n
                diameterInMeters = obj.MICP{i}.PoreThroatDiameter(1:end - 1);
                radiusInMeters = diameterInMeters ./ 2;
                plot(radiusInMeters, slope{i}, 'Linewidth', 2)
                
                [textSlope, textIndex] = max(slope{i});
                text(radiusInMeters(textIndex), textSlope, num2str(obj.MICP{i}.Properties.UserData))                
            end
            xlabel('Pore radius in meters')
            switch stringType
                case 'linear'
                    ylabel('dV/dr')
                case 'log'
                    ylabel('dV/dlog(r)')
            end
            title('Pore Size Distribution')
            axis([1e-9, 1e-4, -inf, inf])
            set(gca, 'XDir', 'reverse')
            set(gca, 'Xscale', 'log')
        end
        
        
        
        function [ pcgwInterp ] = CalcPcgw( obj , nonwettingSaturation )
            pcgwInterp = interp1(obj.MICPInterp.SNW, obj.MICPInterp.PcGW, nonwettingSaturation);
        end
        function [ pchwInterp ] = CalcPchw( obj , nonwettingSaturation )
            pchwInterp = interp1(obj.MICPInterp.SNW, obj.MICPInterp.PcHW, nonwettingSaturation);
        end
        
        function [ averagePcgwPlot , snwPlot ] = GetInterpMICPForPlotting( obj )
            snwPlot = (0:0.01:1)';
            n = numel(snwPlot);
            averagePcgwPlot = zeros(n, 1);
            
            for i = 1:n
                averagePcgwPlot(i) = obj.CalcPcgw(snwPlot(i));
            end
        end
        
        function [ solFigure ] = PlotSol( obj , solFigure , exportTable , doPlotBulkAndMinSol )
            solFigure = obj.PlotSol@BCFormation( solFigure , exportTable , doPlotBulkAndMinSol );
            
            figure(solFigure)
%             axis(obj.solAxis)
            axis(obj.solAxisTemp)
            
            title('Theoretical Formation - Solubility Path')
        end
        function [ sat2PFigure ] = PlotSat2P( obj , sat2PFigure , exportTable , transitionZoneProperties , lineStyle )
            sat2PFigure = obj.PlotSat2P@BCFormation( sat2PFigure , exportTable , transitionZoneProperties , lineStyle );
            
            figure(sat2PFigure)
            axis(obj.satAxis)
            title('Theoretical Formation - 2 Phase Case')
        end
        function [ sat3PFigure ] = PlotSat3P( obj , sat3PFigure , exportTable , lineStyle )
            sat3PFigure = obj.PlotSat3P@BCFormation( sat3PFigure , exportTable , lineStyle );
            
            figure(sat3PFigure)
            axis(obj.satAxis)
            title('Theoretical Formation - 3 Phase Case')
        end        
        function [ pcgwFigure ] = PlotPcgw( obj , pcgwFigure , exportTable , transitionZoneProperties , lineStylePc )
            pcgwFigure = PlotPcgw@BCFormation( obj , pcgwFigure , exportTable , transitionZoneProperties , lineStylePc );
            
            figure(pcgwFigure)

            axis(obj.pcgwAxis)
            title('Theoretical Formation - Gas Overpressure')
        end
        function [ ratioFigure ] = PlotRatio( obj , ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio )
            [ ratioFigure ] = PlotRatio@BCFormation( obj , ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio );
            
            figure(ratioFigure)
            
            title('Theoretical Formation - Overpressure Ratio')
            axis(obj.ratioAxis)
        end
    end
    methods (Static)
        %%% Main function for running 
        function [ minQuantityToFracture3PList , minQuantityToFracture2PList , sgFracture3PList , sgFracture2PList , errorList , depthStructList ] = RunMethaneQuantityFractureRoutine( seafloorDepthArray , initialMethaneQuantity )
            n = numel(seafloorDepthArray);
            errorList = cell(n, 1);
            minQuantityToFracture3PList = zeros(n, 1);
            minQuantityToFracture2PList = zeros(n, 1);
            sgFracture3PList = zeros(n, 1);
            sgFracture2PList = zeros(n, 1);
            depthStructList = cell(n, 1);
            
            for i = 1:n
                iError = 1;
                errorLog = cell(1, 2);
                errorLog{1} = 'SolubilitySaturationRoutine ran and exited without errors';
                found3PFracture = false;
                found2PFracture = false;
                firstTimePlottingSol = true;
                
                transitionZoneProperties.Top3PIndex = [];
                transitionZoneProperties.Bottom3PIndex = [];
                
                seafloorDepth = seafloorDepthArray(i);
                obj = DCTheoreticalFormation(seafloorDepth, 0.4);
                
                if i == 1
                    ch4Quantity = initialMethaneQuantity;
                else
                    ch4Quantity = minQuantityToFracture2PList(i - 1);
                end
                
                
                depthStruct = struct();
                depthStruct.Top3P = [];
                depthStruct.Bottom3P = [];
                depthStruct.BulkEQL3Pfor2P = [];
                
%                 solFigure = figure();
                
                while true
                    try
                        [exportTable, transitionZoneProperties] = obj.RunSolubilitySaturationRoutine(ch4Quantity);
                        exportTable = obj.RunRockAndRatioRoutine(exportTable);
                        
                        if isempty(transitionZoneProperties.Bottom3PIndex)
                            error('Cannot find bottom of three-phase zone')
                        end
                        
                        [sg3PAtFracture, sg2PAtFracture] = obj.CheckFractureStatus(exportTable, transitionZoneProperties);
                        
                        if ~found3PFracture && ~isempty(sg3PAtFracture)
                            minQuantityToFracture3PList(i) = ch4Quantity;
                            sgFracture3PList(i) = sg3PAtFracture;
                            depthStruct.Top3P = exportTable.Depth(transitionZoneProperties.Top3PIndex);
                            depthStruct.Bottom3P = exportTable.Depth(transitionZoneProperties.Bottom3PIndex);
                            found3PFracture = true;
                        end
                        if  ~found2PFracture && ~isempty(sg2PAtFracture)
                            minQuantityToFracture2PList(i) = ch4Quantity;
                            sgFracture2PList(i) = sg2PAtFracture;
                            depthStruct.BulkEQL3Pfor2P = transitionZoneProperties.Bulk3PSolEQLIndex;
                            found2PFracture = true;
                        end
                        
                        
                        
%                         if firstTimePlottingSol
%                             firstTimePlottingSol = false;
%                             doPlotBulkAndMinSol = true;
%                         else
%                             doPlotBulkAndMinSol = false;
%                         end
%                         obj.PlotSolLooped(solFigure, exportTable, transitionZoneProperties, doPlotBulkAndMinSol);
%                         pause(1e-3)
                        
                        
                        
                        
                        input.seafloorDepth = seafloorDepth;
                        input.ch4Quantity = ch4Quantity;
                        input.found3PFracture = found3PFracture;
                        input.found2PFracture = found2PFracture;
                        input.iTop3P = transitionZoneProperties.Top3PIndex;
                        input.iBottom3P = transitionZoneProperties.Bottom3PIndex;
                        DCTheoreticalFormation.PrintRunStatus('RunSolubilitySaturationRoutine', 'Completed successfully', input);
                        
                        
                        
                        if found3PFracture && found2PFracture
                            disp('--------------------------------------------------------')
                            disp('Found both methane quantities to fracture 2P and 3P, breaking search')
                            break
                        end
                    catch exception
                        
                        input.seafloorDepth = seafloorDepth;
                        input.ch4Quantity = ch4Quantity;
                        input.found3PFracture = found3PFracture;
                        input.found2PFracture = found2PFracture;
                        input.iTop3P = transitionZoneProperties.Top3PIndex;
                        input.iBottom3P = transitionZoneProperties.Bottom3PIndex;
                        DCTheoreticalFormation.PrintRunStatus('RunSolubilitySaturationRoutine', 'Exited on error', input);
                        
                        [errorLog, iError, executeBreak] = DCTheoreticalFormation.RunErrorRoutine(exception, ch4Quantity, errorLog, iError);
                        
                        if executeBreak
                            break
                        end
                        
                    end
                    
                    ch4Quantity = ch4Quantity + 1;
                end
                
                depthStructList{i} = depthStruct;
                errorList{i} = errorLog;
            end
        end
        
        %%% Plotting functions
        function PlotMethaneQuantities( seafloorDepthArray , minQuantityToFracture2PList , minQuantityToFracture3PList )
            figure
            hold on
            plot(seafloorDepthArray, minQuantityToFracture2PList, 'linewidth', 3)
            plot(seafloorDepthArray, minQuantityToFracture3PList, 'linewidth', 3)
            xlabel('Seafloor depth (m)')
            ylabel('Minimum methane quantity to initiate a fracture (kg/m^3)')
            legend('2P case', '3P case')
            
            figure
            hold on
            plot(seafloorDepthArray, minQuantityToFracture3PList ./ minQuantityToFracture2PList, 'linewidth', 3)
            xlabel('Seafloor depth (m)')
            ylabel('Ratio of 3P over 2P required methane quantity to initiate a fracture')
        end
        function PlotGasSaturations( seafloorDepthArray , sgFracture2PList , sgFracture3PList )
            figure
            hold on
            plot(seafloorDepthArray, sgFracture2PList, 'linewidth', 3)
            plot(seafloorDepthArray, sgFracture3PList, 'linewidth', 3)
            xlabel('Seafloor depth (m)')
            ylabel('Gas Saturation at fracture initiation')
            legend('2P case', '3P case')
        end
        function PlotDepths( seafloorDepthArray , depthStructList )
            n = numel(depthStructList);
            
            top3PDepths = nan(n, 1);
            bottom3PDepths = nan(n, 1);
            bulkEQL3PDepths = nan(n, 1);
            
            for i = 1:n
                depthStruct = depthStructList{i};
                if ~isempty(depthStruct.BulkEQL3Pfor2P)
                    bulkEQL3PDepths(i) = depthStruct.BulkEQL3Pfor2P;
                end
                
                if ~isempty(depthStruct.Top3P)
                    top3PDepths(i) = depthStruct.Top3P;
                end
                if ~isempty(depthStruct.Bottom3P)
                    bottom3PDepths(i) = depthStruct.Bottom3P;
                end
            end
            
            
            figure
            hold on
%             scatter(seafloorDepthArray, bulkEQL3PDepths, 50)
            scatter(seafloorDepthArray, bulkEQL3PDepths ./ 2, 50) % hotfix because i didn't save the depths @@@ NEED TO FIX
            plot(seafloorDepthArray, top3PDepths, 'linewidth', 3)
            plot(seafloorDepthArray, bottom3PDepths, 'linewidth', 3)
            xlabel('Seafloor depth (m)')
            ylabel('Depth (mbsf)')
            legend('2P Case - Bulk equilibrium depth at fracture', '3P case - Top of 3P zone at fracture', '3P case - Bottom of 3P zone at fracture')
            set(gca, 'YDir', 'Reverse')
        end        
        
        
        %%% Error and logging functions
        function PrintRunStatus( functionString , message , input )
            
            disp('--------------------------------------------------------')
            disp(['Running: ', functionString])
            disp(['Current seafloor depth (m): ', num2str(input.seafloorDepth)])
            disp(['Current methane quantity (kg/m^3): ', num2str(input.ch4Quantity)])
            disp(['3P fracture found: ', char(string(input.found3PFracture))])
            disp(['2P fracture found: ', char(string(input.found2PFracture))])
            if isempty(input.iTop3P)
                disp('Index of 3P zone top: not found')
            else
                disp(['Index of 3P zone top: ', num2str(input.iTop3P)])
            end
            if isempty(input.iBottom3P)
                disp('Index of 3P zone bottom: not found')
            else
                disp(['Index of 3P zone bottom: ', num2str(input.iBottom3P)])
            end
            disp(message)
            
        end
        function [ errorLog , iError , executeBreak ] = RunErrorRoutine( exception , ch4Quantity , errorLog , iError )
            executeBreak = false;
            replaceMessageWithException = false;
            
            isSg = contains(lower(exception.message), 'maxsollg');
            isSh = contains(lower(exception.message), 'maxsollh');
            is3P = contains(lower(exception.message), '3p');
            
            startIndex = regexp(exception.message, '[0-9]\.');
            if isempty(startIndex)
                errorNumber = -1;
                string2 = string();
            else
                errorNumber = str2double(exception.message(startIndex:end));
                string2 = num2str(errorNumber);
            end
            
            if isSg
                string1 = 'Sg is: ';
            elseif isSh
                string1 = 'Sh is: ';
            elseif is3P
                string1 = '3P calculation error';
            else
                string1 = '';
                replaceMessageWithException = true;
            end
            
            errorLog{iError, 2} = ch4Quantity;
            errorLog{iError, 1} = [string1, string2];
            if replaceMessageWithException
                errorLog{iError, 1} = exception;
            end
            
            if errorNumber > 1
                disp('--------------------------------------------------------')
                disp('Found max methane quantity (2P saturation > 1), breaking search')
                executeBreak = true;
            end
            iError = iError + 1;
        end
        
        
        
        
        %%% MICP functions
        function ExtractMICP()
            % Used to assign each MICP report with the core's drilled depth
            depth =  {  'TABLE_S1_KZK867.XLSX',	925.5;
                        'TABLE_S2_ERIS89.XLSX',	1320.5;
                        'TABLE_S3_3NIMU5.XLSX',	1555.5;
                        'TABLE_S4_4XF4BE.XLSX',	1725.5;
                        'TABLE_S5_WVC65Q.XLSX',	1925.5;
                        'TABLE_S6_BQQWYW.XLSX',	1977.5;
                        'TABLE_S7_4DNM5V.XLSX',	1110.5;
                        'TABLE_S8_V9IEVD.XLSX',	905;
                        'TABLE_S9_J0BZ4M.XLSX',	912;
                        'TABLE_S10_D1IQ8N.XLSX',	925;
                        'TABLE_S11_Z0S52K.XLSX',	221.25;
                        'TABLE_S12_X8P4Y3.XLSX',	253.75;
                        'TABLE_S13_42778W.XLSX',	280;
                        'TABLE_S14_DXWGHA.XLSX',	311.5;
                        'TABLE_S15_ACU9AK.XLSX',	348.75;
                        'TABLE_S16_DAL6FO.XLSX',	376;
                        'TABLE_S17_RDWGHU.XLSX',	405.25;
                        'TABLE_S18_KLXRT5.XLSX',	435;
                        'TABLE_S19_O73KNK.XLSX',	463.5   };
            
            
            
            % Change directories and get file names
            oldDir =  cd('Kumano Basin Data\Mercury Intrusion C0002\');
            allFiles = dir;
            allNames = { allFiles(~[allFiles.isdir]).name }';
            allNamesLogical = contains(allNames, '.XLSX');
            allNames = allNames(allNamesLogical);
            nFile = numel(allNames);
            if nFile ~= 19
                error('Number of Excel files does not match 19')
            end
            
            % Prep MICP data gather
            MICPCellArray = cell(nFile, 1);
            
            sigmaHgAir = 485; % dynes/cm
            sigmaGasWater = 72; % dynes/cm
            sigmaHydrateWater = 27; % dynes/cm
            
            thetaHgAir = 140; % deg
            thetaGasWater = 180; % deg
            thetaHydrateWater = 180; % deg
            
            conversionPsiToMpa = 1/145.03774; % multiply by this
            
            % Forloop through all xlsx in data directory
            for iFile = 1:nFile
                fileString = allNames{iFile};
                tempMatrix = xlsread(fileString);
                
                if size(tempMatrix, 2) ~= 4
                    disp(allNames{iFile});
                    error('Excel file does not have 4 columns of data');
                end
                MICP = table();
                
                % Save non-wetting saturation
                MICP.SNW = tempMatrix(:, 2);
                
                % Save Pcgw and Pchw converting from Pc mercury-air
                MICP.PcGW = tempMatrix(:, 1) .* (sigmaGasWater * cosd(thetaGasWater)) ...
                                              / (sigmaHgAir * cosd(thetaHgAir)); % psi
                MICP.PcHW = tempMatrix(:, 1) .* (sigmaHydrateWater * cosd(thetaHydrateWater)) ...
                                              / (sigmaHgAir * cosd(thetaHgAir)); % psi
                MICP.PcGW = MICP.PcGW .* conversionPsiToMpa;
                MICP.PcHW = MICP.PcHW .* conversionPsiToMpa;
                
                % Convert pore throat radius in microns to pore throat diameter in meters
                MICP.PoreThroatDiameter = tempMatrix(:, 4) .* 2 ./ 1e6;
                
                
                fileLogical = contains(depth(:, 1), fileString);
                MICP.Properties.UserData = depth{fileLogical, 2};
                
                
                MICPCellArray{iFile} = MICP;
            end
            cd(oldDir)
            save('MICP_KB.mat', 'MICPCellArray')
        end
        function [ result ] = LoadMICP()
            MICPCellArray = [];
            % Loads MICPCellArray
            load('MICP_KB.mat');
%             
%             n = numel(MICPCellArray);
%             depthOrder = zeros(n, 1);
%             for i = 1:n
%                 depthOrder(i) = MICPCellArray{i}.Properties.UserData;
%             end
%             [~, colorOrder] = sort(depthOrder);
%             result = MICPCellArray{8}.Properties.UserData;
            result = MICPCellArray{8};
        end
        function [ MICPInterp ] = SelectMICPForInterp( MICP )
%             logicalToKeep = depthOrder > 340 & depthOrder < 500;
%             logicalToKeep = depthOrder > 400 & depthOrder < 420;
            MICPInterp = MICP;
            
%             MICPInterp(~logicalToKeep) = [];
%             n = numel(MICPInterp);
            
%             for i = 1:n
            temp = MICPInterp;
            firstNonZeroIndex = find(temp.SNW > 0, 1);
            temp(1:firstNonZeroIndex - 2, :) = [];

            MICPInterp = temp;
%             end
            
        end
        
    end
    % UNUSED CLASS METHODS
    %{
        function [ bulkDensity , porosity ] = EstimateBulkDensity1( obj )
            % Including the non-logged depths (in mbsf) in the effective vertical stress
            % This function is only used for the fracture code
            
            depth = obj.DataTable.depth;
            
            
            Phi_0 = obj.phi0;
            Phi_inf = obj.phiInf;
            B = obj.depthB; % meters
            
            Rho_fluid = 1.024; % g/cc, seawater
            Rho_grain = 2.7;   % g/cc, smectite
            
            porosity = Phi_inf + (Phi_0 - Phi_inf)*exp(-depth./B);
            bulkDensity = porosity*Rho_fluid + (1 - porosity)*Rho_grain;
        end
        function [ ratioFigure ] = PlotFractureRatio1( ~ , ratioFigure )
            figure(ratioFigure)
            
            plot( [1 1] , [0 300] , 'k--' , 'linewidth' , 2 )
            hold on
        end
        function LoadSolubility( obj )
            load('Blake Ridge Data\Solubility Plots\HR_bulk_C_L_G.mat')
            obj.Bulk.Solubility = HR_bulk_C_L_G(:,1); % mol CH4/kg H2O
            obj.Bulk.Depth = HR_bulk_C_L_G(:,2) - obj.seafloorDepth; % mbsf
                        
        end
        function PreprocessData( obj )
            
            obj.Data.depth_top = obj.Depth.saturationTop;
            obj.Data.depth_bottom = obj.Depth.saturationBottom;
            
            obj.Saturation.Above_Data = repmat( obj.Saturation.HydrateGrid(1,:) , obj.Depth.Graph_Top - obj.Depth.saturationTop , 1 );
            obj.Saturation.Below_Data = repmat( obj.Saturation.GasGrid(end,:) , obj.Depth.saturationBottom - obj.Depth.Graph_Bottom , 1 );
            
            obj.Saturation.Hydrate_Data = vertcat(obj.Saturation.Above_Data,obj.Saturation.HydrateGrid);
            obj.Saturation.Gas_Data = vertcat(obj.Saturation.GasGrid,obj.Saturation.Below_Data);
            

            obj.Data.depth_interval = obj.Data.depth_bottom - obj.Data.depth_top + 1;
            obj.Data.log = zeros( obj.Data.depth_interval , 50 );
            obj.Data.log(:,1) = linspace( obj.Data.depth_top , obj.Data.depth_bottom , obj.Data.depth_interval );

%             obj.DataTable.depth = (obj.Data.depth_top : 1 : obj.Data.depth_bottom)';
            obj.DataTable.depth = (obj.Data.depth_top : 0.5 : obj.Data.depth_bottom)';

        end
    %}
    % UNUSED STATIC METHODS
    %{
        function [ result ] = LoadPhaseBehaviorHydrateRidge()
            % Loads saturation data from hydrate and gas grabit array .mat files and
            %   saves into struct Saturation with fields 'Hydrate' and 'Gas'
            Saturation = struct('Hydrate',cell2mat(struct2cell(Load('Hydrate Ridge Data\Phase Behavior\HydrateBehavior.mat'))),'Gas',cell2mat(struct2cell(Load('Hydrate Ridge Data\Phase Behavior\GasBehavior.mat'))));
            
            % Upper and lower depths of grabit depth of investigation
            Depth_grid = linspace(880,950,71);
            Methane_quantity_grid = linspace(0,40,41);
            
            [X,Y] = meshgrid(Methane_quantity_grid,Depth_grid);
            
            Saturation.HydrateGrid = griddata(Saturation.Hydrate(:,1),Saturation.Hydrate(:,2),Saturation.Hydrate(:,3),X,Y,'cubic');
            Saturation.HydrateGrid(Saturation.HydrateGrid < 0) = 0;
            Saturation.HydrateGrid(isnan(Saturation.HydrateGrid)) = 0;
            Saturation.HydrateGrid(Saturation.HydrateGrid < .0001 & Saturation.HydrateGrid ~= 0) = 0;
            
            Saturation.GasGrid = griddata(Saturation.Gas(:,1),Saturation.Gas(:,2),Saturation.Gas(:,3),X,Y,'cubic');
            Saturation.GasGrid(Saturation.GasGrid < 0) = 0;
            Saturation.GasGrid(isnan(Saturation.GasGrid)) = 0;
            Saturation.GasGrid(Saturation.GasGrid < .0001 & Saturation.GasGrid ~= 0) = 0;
            
            result = Saturation;
        end
    
        function FitBulkDensityData()
            %%% NOT FINISHED, NOT USED
            %%% use the excel in \Kumano Basin Data\
            
            KumanoBasinData = [];
            % Loads KumanoBasinData
            load('RHOB_KB.mat');
            depth = KumanoBasinData.Depth;
            bulkDensity = KumanoBasinData.RHOB;
            
            
            
            
            
            expFitType = fittype('a + b*log(x)', ...
                'dependent', {'y'}, ...
                'independent', {'x'}, ...
                'coefficients', {'a', 'b'});
            [logFit logGOF] = fit(xData, yData, expFitType)
            
            
            
            figure
            hold on
            
            xlabel('Bulk density')
            ylabel('Depth (mbsf)')
            set(gca, 'YDir', 'Reverse')
            
            
        end

    %}
end























