classdef DCKumanoBasin < BCFormation
    properties
        MICP % Cell array of tables
        colorOrder % array containing order index of loaded MICP tables
        depthOrder % array containing depths in the order of loaded MICP tables
        MICPInterp % Cell array of tables of MICP that will actually be used
    end
    properties (Constant)
        satAxis = [0 1 400 430];
        solAxis = [0.155 0.17 400 430];
        pcgwAxis = [0 2.5 400 430];
        ratioAxis = [0 1.2 400 430];
        

    end
    methods
        %%% Constructor
        function [ obj ] = DCKumanoBasin()
            obj@BCFormation()
            
            obj.seafloorDepth = 1936;        % m
            obj.minDepth = 1;               % mbsf
            obj.maxDepth = 500;             % mbsf
            obj.depthIncrement = 0.5;       % m
            obj.depthArray = (obj.minDepth : obj.depthIncrement : obj.maxDepth)';
            
            obj.temperatureGradient = 40;   % C deg/km (C0002 - Daigle and Dugan 2014)            
            obj.seafloorTemperature = 2.2;    % C deg
            obj.salinityWtPercent = 3.5;    % weight percent (wt%) of NaCl in seawater
            
            obj.phi0 = 0.6;
            obj.phiInf = 0.25;
            obj.B = 1000; % m
            
            % Loads MICP from a .mat file into the object property MICP
            [obj.MICP, obj.colorOrder, obj.depthOrder] = DCKumanoBasin.LoadMICP();
            obj.MICPInterp = DCKumanoBasin.SelectMICPForInterp(obj.MICP, obj.depthOrder);
            
            % Testing plots only on selected MICP data for interpolation
            obj.MICP = obj.MICPInterp;
        end
        
        %%% Petrophysical calculations
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
        function [ pcgwInterp ] = CalcPcgw( obj , nonwettingSaturation )
            % accepts nonwettingSaturation as an array
            nInterpSets = numel(obj.MICPInterp);
            nSnw = numel(nonwettingSaturation);
            
            pcgwSets = zeros(nSnw, nInterpSets);
            
            
            for iInterpSets = 1:nInterpSets
                tempTable = obj.MICPInterp{iInterpSets};
                pcgwSets(:, iInterpSets) = interp1(tempTable.SNW, tempTable.PcGW, nonwettingSaturation);
            end
            pcgwInterp = mean(pcgwSets, 2);
        end
        function [ pchwInterp ] = CalcPchw( obj , nonwettingSaturation )
            % does not accept nonwettingSaturation as an array
            n = numel(obj.MICPInterp);
            pcgwArray = zeros(n, 1);
            for i = 1:n
                tempTable = obj.MICPInterp{i};
                pcgwArray(i) = interp1(tempTable.SNW, tempTable.PcHW, nonwettingSaturation);
            end
            pchwInterp = mean(pcgwArray);
        end
        
        function [ averagePcgwPlot , snwPlot ] = GetInterpMICPForPlotting( obj )
            snwPlot = (0:0.01:1)';
            n = numel(snwPlot);
            averagePcgwPlot = zeros(n, 1);
            
            for i = 1:n
                averagePcgwPlot(i) = obj.CalcPcgw(snwPlot(i));
            end
        end
        
        function [ solFigure ] = PlotSol( obj , solFigure , exportTable )
            solFigure = obj.PlotSol@BCFormation( solFigure , exportTable );
            
            figure(solFigure)
            axis(obj.solAxis)
            title('Kumano Basin - Solubility Path')
        end
        function [ sat2PFigure ] = PlotSat2P( obj , sat2PFigure , exportTable , transitionZoneProperties , lineStyle )
            sat2PFigure = obj.PlotSat2P@BCFormation( sat2PFigure , exportTable , transitionZoneProperties , lineStyle );
            
            figure(sat2PFigure)
            axis(obj.satAxis)
            title('Kumano Basin - 2 Phase Case')
        end
        function [ sat3PFigure ] = PlotSat3P( obj , sat3PFigure , exportTable , lineStyle )
            sat3PFigure = obj.PlotSat3P@BCFormation( sat3PFigure , exportTable , lineStyle );
            
            figure(sat3PFigure)
            axis(obj.satAxis)
            title('Kumano Basin - 3 Phase Case')
        end        
        function [ pcgwFigure ] = PlotPcgw( obj , pcgwFigure , exportTable , transitionZoneProperties , lineStylePc )
            pcgwFigure = PlotPcgw@BCFormation( obj , pcgwFigure , exportTable , transitionZoneProperties , lineStylePc );
            
            figure(pcgwFigure)

            axis(obj.pcgwAxis)
            title('Kumano Basin Gas Overpressure')
            % legend('')
        end
        function [ ratioFigure ] = PlotRatio( obj , ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio )
            [ ratioFigure ] = PlotRatio@BCFormation( obj , ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio );
            
            figure(ratioFigure)
            
            title('Hydrate Ridge Overpressure Ratio')
            axis(obj.ratioAxis)            
            % legend('')
            
        end
    end
    methods (Static)
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
        function [ result , colorOrder , depthOrder ] = LoadMICP()
            MICPCellArray = [];
            % Loads MICPCellArray
            load('MICP_KB.mat');
            
            n = numel(MICPCellArray);
            depthOrder = zeros(n, 1);
            for i = 1:n
                depthOrder(i) = MICPCellArray{i}.Properties.UserData;
            end
            [~, colorOrder] = sort(depthOrder);
            result = MICPCellArray;
        end
        function [ MICPInterp ] = SelectMICPForInterp( MICP , depthOrder )
%             logicalToKeep = depthOrder > 340 & depthOrder < 500;
            logicalToKeep = depthOrder > 400 & depthOrder < 420;
            MICPInterp = MICP;
            
            MICPInterp(~logicalToKeep) = [];
            n = numel(MICPInterp);
            
            for i = 1:n
                temp = MICPInterp{i};
                firstNonZeroIndex = find(temp.SNW > 0, 1);
                temp(1:firstNonZeroIndex - 2, :) = [];
                
                MICPInterp{i} = temp;
            end
            
        end

    end
    % UNUSED CLASS METHODS
    %{
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























