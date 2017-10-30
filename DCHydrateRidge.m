classdef DCHydrateRidge < BCFormation
    properties
        MICP1
    end
    properties (Constant)
        satAxis = [0 1 100 160];
    end
    methods
        %%% Constructor
        function [ obj ] = DCHydrateRidge()
            obj@BCFormation()
            
            obj.seafloorDepth = 790;        % m
            obj.minDepth = 1;               % mbsf
            obj.maxDepth = 200;             % mbsf
            obj.depthIncrement = 0.5;       % m
            obj.depthArray = (obj.minDepth : obj.depthIncrement : obj.maxDepth)';
            
            obj.temperatureGradient = 59;   % C deg/km (ODP 997 - Liu and Flemings)            
            obj.seafloorTemperature = 4;    % C deg
            obj.salinityWtPercent = 3.5;    % weight percent (wt%) of NaCl in seawater
            
            obj.MICP1 = DCHydrateRidge.LoadMICP1();       
        end
        
        %%% Petrophysical calculations
        function [ bulkDensity , porosity ] = EstimateBulkDensity( obj )
            % Including the non-logged depths (in mbsf) in the effective vertical stress
            % This function is only used for the fracture code
            
            depth = obj.DataTable.depth;
            
            
            Phi_0 = 0.63;
            Phi_inf = 0.1;
            B = 1400; % meters
            
            Rho_fluid = 1.024; % g/cc, seawater
            Rho_grain = 2.7;   % g/cc, smectite
            
            porosity = Phi_inf + (Phi_0 - Phi_inf)*exp(-depth./B);
            bulkDensity = porosity*Rho_fluid + (1 - porosity)*Rho_grain;            
        end
        function [ pcgwInterp ] = CalcPcgw( obj , nonwettingSaturation )
            pcgwInterp = interp1( obj.MICP1.S_nw , obj.MICP1.Pc_gw , nonwettingSaturation );
        end
        function [ pchwInterp ] = CalcPchw( obj , nonwettingSaturation )
            pchwInterp = interp1( obj.MICP1.S_nw , obj.MICP1.Pc_hw , nonwettingSaturation );
        end
        
        function [ diameter ] = CalcPoreDiameterFromPcgw( obj )
            % pcgw in MPa
            % diameter in microns (1e-6 m)
            
            pcgw = obj.MICP1.Pc_gw;
            diameter = 4 * 0.072 ./ pcgw;
        end
        function [ slope ] = CalcSlopeOfCumPSD( obj, stringType )
            % pcgw in MPa
            % diameter in microns (1e-6 m)
            
            pcgw = obj.MICP1.Pc_gw;
            diameter = 4 * 0.072 ./ pcgw;
            sNw = obj.MICP1.S_nw;
            
            switch stringType
                case 'linear'
                    slope =     -( sNw(2:end) - sNw(1:end - 1) ) ...
                              ./ ( diameter(2:end) - diameter(1:end - 1) );
                    
                          
%                     slope = zeros(numel(sNw) - 1, 1);
%                     for i = 1:numel(sNw) - 1
%                         slope(i) = -( sNw(i + 1) - sNw(i) ) ...
%                                   / (diameter(i + 1) - diameter(i) );
%                         
%                         
%                     end
                case 'log'
                    slope =     -( sNw(2:end) - sNw(1:end - 1) ) ...
                              ./ ( log10(diameter(2:end)) - log10(diameter(1:end - 1)) );
            end
        end
        
        
        
        %%% Plotting subclass functions
        function [ solFigure ] = PlotSol( obj , solFigure , exportTable )
            solFigure = obj.PlotSol@BCFormation( solFigure , exportTable );
            
            figure(solFigure)
            axis([0.08 0.125 60 160])
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
            
            [ diameter ] = obj.CalcPoreDiameterFromPcgw();
            % diameter is in microns, plot converts to meters
            figure
            semilogx(diameter ./ 1e6, obj.MICP1.S_nw, 'Linewidth', 2)
%             xlabel('Pore diameter in microns (1e-6 m)')
            xlabel('Pore diameter in meters')
            ylabel('S_n_w')
            title('Cumulative Pore Size Distribution')
            axis([1e-9, 1e-4, 0, 1])
            set(gca, 'XDir', 'reverse')
        end
        function PlotPSD( obj , stringType )
            
            diameter = obj.CalcPoreDiameterFromPcgw();
            slope = obj.CalcSlopeOfCumPSD( stringType );
            
            
            % diameter is in microns, plot converts to meters
            figure
            semilogx(diameter(1:end - 1) ./ 1e6, slope, 'Linewidth', 2)
            xlabel('Pore diameter in meters')
%             ylabel('Frequency f(r)')
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
        
        % work on these next

        function [ pcgwFigure ] = PlotPcgw( obj , pcgwFigure , iStorage , lineStyle2D , lineStyle3D )
            [ pcgwFigure ] = PlotPcgw@Formation( obj , pcgwFigure , iStorage , lineStyle2D , lineStyle3D );
            
            figure(pcgwFigure)

            axis([0 2 100 150])
            title('Hydrate Ridge Gas Overpressure')
            % legend('')
        end
        function [ ratioFigure ] = PlotRatio( obj , ratioFigure , iStorage , lineStyle2D , lineStyle3D )
            [ ratioFigure ] = PlotRatio@Formation( obj , ratioFigure , iStorage , lineStyle2D , lineStyle3D );
            
            figure(ratioFigure)
            
            title('Hydrate Ridge Overpressure Ratio')
            axis([0 2.5 100 150])            
            % legend('')
            
        end
        function [ ratioFigure ] = PlotFractureRatio( ~ , ratioFigure )
            figure(ratioFigure)
            
            plot( [1 1] , [0 300] , 'k--' , 'linewidth' , 2 )
            hold on
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

    %}
end















































