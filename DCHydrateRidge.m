classdef DCHydrateRidge < BCFormation
    properties
        MICP1
    end
    properties (Constant)
        poissonsRatio = 0.4;
        
        satAxis = [0 1 100 160];
        solAxis = [0.1 0.125 100 160];
        pcgwAxis = [0 2 100 160];
        ratioAxis = [0 3.5 100 160];
    end
    methods
        %%% Constructor
        function [ obj ] = DCHydrateRidge()
            obj@BCFormation()
            
            obj.seafloorDepth = 790;        % m
            obj.minDepth = 1;               % mbsf
            obj.maxDepth = 200;             % mbsf
            obj.depthIncrement = 0.5;       % m
            % obj.depthIncrement = 0.25;       % m
            obj.depthArray = (obj.minDepth : obj.depthIncrement : obj.maxDepth)';
            
            obj.temperatureGradient = 59;   % C deg/km (ODP 997 - Liu and Flemings)            
            obj.seafloorTemperature = 4;    % C deg
            obj.salinityWtPercent = 3.5;    % weight percent (wt%) of NaCl in seawater
            
            obj.phi0 = 0.63;
            obj.phiInf = 0.1;
            obj.B = 1400; % m
            
            
            
            obj.MICP1 = DCHydrateRidge.LoadMICP1();
            [obj.sNwArray, obj.radiusArray, obj.lengthArray, obj.nArray] = obj.CalcPoreVolumeDistribution();
        end
        
        %%% Petrophysical calculations
        function [ pcgwInterp ] = CalcPcgw( obj , nonwettingSaturation )
            pcgwInterp = interp1( obj.MICP1.S_nw , obj.MICP1.Pc_gw , nonwettingSaturation );
        end
        function [ pchwInterp ] = CalcPchw( obj , nonwettingSaturation )
            pchwInterp = interp1( obj.MICP1.S_nw , obj.MICP1.Pc_hw , nonwettingSaturation );
        end
        
        function [ radiusInMeters ] = CalcPoreRadiusFromPcgw( obj )
            % pcgw in MPa
            % radius in meters
            
            pcgw = obj.MICP1.Pc_gw;
            radiusInMicrons = 2 * 0.072 ./ pcgw;
            radiusInMeters = radiusInMicrons ./ 1e6;
        end
        function [ slope ] = CalcSlopeOfCumPSD( obj, stringType )
            % slope is in volume fraction per radius meter
            
            radiusInMeters = obj.CalcPoreRadiusFromPcgw();
            sNw = obj.MICP1.S_nw;
            
            switch stringType
                case 'linear'
                    slope =     -( sNw(2:end) - sNw(1:end - 1) ) ...
                              ./ ( radiusInMeters(2:end) - radiusInMeters(1:end - 1) );
                case 'log'
                    slope =     -( sNw(2:end) - sNw(1:end - 1) ) ...
                              ./ ( log10(radiusInMeters(2:end)) - log10(radiusInMeters(1:end - 1)) );
            end
        end
        function [ sg , radiusGPrime ] = InterpCumPSD( obj , shGuess , radiusG , radiusH )
            radiusInMetersArray = obj.CalcPoreRadiusFromPcgw();
            snwArray = obj.MICP1.S_nw;
            
            snw = interp1(radiusInMetersArray, snwArray, radiusG);
            sw = 1 - snw;
            
            sg1 = snw - interp1(radiusInMetersArray, snwArray, radiusH);
            
            sg2 = 1 - sw - shGuess - sg1;
            
            sg = sg1 + sg2;
            radiusGPrime = interp1(snwArray, radiusInMetersArray, sg2);
            
            % scatter(radiusG, snw, 'r', 'filled')
            % scatter(radiusH, shGuess + sg2, 'g', 'filled')
            % scatter(radiusGPrime, sg2, 'r', 'filled')
            
            
            
            
        end
        
        function [ sNwArray , radiusArray , lengthArray , nArray ] = CalcPoreVolumeDistribution( obj )
            

            numberOfPoints = 1000;
            
            slopeScalingFactor = (1e-7) ^ 2;
            lengthFactor = 3;
            
            radiusInMeters = obj.CalcPoreRadiusFromPcgw();
            
            %%% Units are in meters
            radiusMin = min(log10(radiusInMeters));
            radiusMax = max(log10(radiusInMeters));
            radiusArray = logspace(radiusMax, radiusMin, numberOfPoints)';
            
            sNwArray = interp1(radiusInMeters, obj.MICP1.S_nw, radiusArray);
            if isnan(sNwArray(end))
                sNwArray(end) = 1;
            end
            dVdrArray =   -( sNwArray(2:end) - sNwArray(1:end - 1) ) ...
                        ./ ( radiusArray(2:end) - radiusArray(1:end - 1) );
            % makes slope array same length as the rest of the arrays
            % also scales the slope to something reasonable instead of
            % assuming the actual volume change is on the order of 1 m^3
            dVdrArray = [dVdrArray; dVdrArray(end)] .* slopeScalingFactor;
            
            
            lengthArray = radiusArray .* lengthFactor;
            
            % derived expression assuming lengthFactor = 3
            nArray = dVdrArray ./ (9 .* pi .* radiusArray .^ 2);
            
            
            
            
            
            return
            
            figure
            semilogx(radiusArray, dVdrArray, 'Linewidth', 2)
            xlabel('Pore radius in meters')
            ylabel('dV/dr')
            axis([1e-9, 1e-4, -inf, inf])
            set(gca, 'XDir', 'reverse')

            figure
            loglog(radiusArray, nArray, 'Linewidth', 2)
            xlabel('Pore radius in meters')
            ylabel('n')
            axis([1e-9, 1e-4, -inf, inf])
            set(gca, 'XDir', 'reverse')
        end
        
        
        %%% Plotting subclass functions
        function [ solFigure ] = PlotSol( obj , solFigure , exportTable , doPlotBulkAndMinSol )
            solFigure = obj.PlotSol@BCFormation( solFigure , exportTable , doPlotBulkAndMinSol );
            
            figure(solFigure)
            axis(obj.solAxis)
            title('Hydrate Ridge - Solubility Path')
        end
        function [ sat2PFigure ] = PlotSat2P( obj , sat2PFigure , exportTable , transitionZoneProperties , lineStyle )
            sat2PFigure = obj.PlotSat2P@BCFormation( sat2PFigure , exportTable , transitionZoneProperties , lineStyle );
            
            figure(sat2PFigure)
            axis(obj.satAxis)
            title('Hydrate Ridge - 2 Phase Case')
        end
        function [ sat3PFigure ] = PlotSat3P( obj , sat3PFigure , exportTable , lineStyle )
            sat3PFigure = obj.PlotSat3P@BCFormation( sat3PFigure , exportTable , lineStyle );
            
            figure(sat3PFigure)
            axis(obj.satAxis)
            title('Hydrate Ridge - 3 Phase Case')
        end        
        
        function PlotMICP( obj )
            figure
            semilogy(1 - obj.MICP1.S_nw, obj.MICP1.Pc_gw, 'Linewidth', 2)
            xlabel('1 - S_n_w, or S_w')
            ylabel('Pc in MPa')
            title('Primary Drainage Capillary Pressure Curve')
        end
        function PlotCumPSD( obj )
            
            radiusInMeters = obj.CalcPoreRadiusFromPcgw();
            
            
            figure
            semilogx(radiusInMeters, obj.MICP1.S_nw, 'Linewidth', 2)
            xlabel('Pore radius in meters')
            ylabel('S_n_w')
            title('Cumulative Pore Size Distribution')
            axis([1e-9, 1e-4, 0, 1])
            set(gca, 'XDir', 'reverse')
        end
        function PlotPSD( obj , stringType )
            
            radiusInMeters = obj.CalcPoreRadiusFromPcgw();
            slope = obj.CalcSlopeOfCumPSD(stringType);
            
            
            figure
            semilogx(radiusInMeters(1:end - 1), slope, 'Linewidth', 2)
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
        end
        function pcgwFigure = PlotPcgw( obj , pcgwFigure , exportTable , transitionZoneProperties , lineStylePc )
            pcgwFigure = PlotPcgw@BCFormation( obj , pcgwFigure , exportTable , transitionZoneProperties , lineStylePc );
            
            figure(pcgwFigure)

            axis(obj.pcgwAxis)
            title('Hydrate Ridge - Gas Overpressure')
        end
        function ratioFigure = PlotRatio( obj , ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio )
            [ ratioFigure ] = PlotRatio@BCFormation( obj , ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio );
            
            figure(ratioFigure)
            
            title('Hydrate Ridge - Overpressure Ratio')
            axis(obj.ratioAxis)
        end
        
    end
    methods (Static)
        function [ result ] = LoadMICP1()
            MICP_HR.file_ID = fopen('MICP_HR.dat');
            MICP_HR.ScannedData = fscanf(MICP_HR.file_ID,'%f',[3 Inf])';
            fclose(MICP_HR.file_ID);
            
            MICP_HR.S_nw = MICP_HR.ScannedData(:,1);
            MICP_HR.Pc_gw = MICP_HR.ScannedData(:,2);
            MICP_HR.Pc_hw = MICP_HR.ScannedData(:,3);
            
            result = MICP_HR;
        end
    end
    % UNUSED CLASS METHODS
    %{
        function [ bulkDensity , porosity ] = EstimateBulkDensity1( obj )
            % Including the non-logged depths (in mbsf) in the effective vertical stress
            % This function is only used for the fracture code
            depth = obj.depthArray;
            
            
            Phi_0 = 0.63;
            Phi_inf = 0.1;
            B = 1400; % meters
            
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

             obj.DataTable.depth = (obj.Data.depth_top : 1 : obj.Data.depth_bottom)';
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

    %}
end















































