classdef DCBlakeRidge < BCFormation
    properties
        LDT
        DIT
        GR
        BRG
        CALI
        
        MICP1
        MICP2
    end
    properties (Constant)
        satAxis = [0 1 460 510];
        solAxis = [0.18 0.19 460 510];
        pcgwAxis = [0 2 460 510];
        ratioAxis = [0 1.2 460 510];
    end
    methods
        %%% Constructor
        function [ obj ] = DCBlakeRidge()
            obj@BCFormation()

            obj.seafloorDepth = 2780;           % mbsl
            obj.minDepth = 1;                   % mbsf
%             obj.maxDepth = 650;                 % mbsf
            obj.maxDepth = 700;                 % mbsf
            obj.depthIncrement = 0.5;           % m
            obj.depthArray = (obj.minDepth : obj.depthIncrement : obj.maxDepth)';
            
%             obj.temperatureGradient = 38.5;     % C deg/km (ODP 995 - Guerin)
            obj.temperatureGradient = 36.9;     % C deg/km (ODP 997 - Liu and Flemings)
            obj.seafloorTemperature = 3.3;      % C deg
            obj.salinityWtPercent = 3.5;        % weight percent (wt%) of NaCl in seawater
            
            obj.phi0 = 0.75;
            obj.phiInf = 0.05;
            obj.B = 1600; % m
            
            
            obj.LDT = DCBlakeRidge.LoadLDT();
            obj.DIT = DCBlakeRidge.LoadDIT();
            obj.GR = DCBlakeRidge.LoadGR();
            obj.BRG = DCBlakeRidge.LoadBRG();
            obj.CALI = DCBlakeRidge.LoadCALIHLDT();
            
            obj.MICP1 = DCBlakeRidge.LoadMICP1();
            obj.MICP2 = DCBlakeRidge.LoadMICP2();
        end
        
        %%% Petrophysical calculations
        function [ bulkDensity , porosity ] = EstimateBulkDensity1( obj )
            % Including the non-logged depths (in mbsf) in the effective vertical stress            
            depth = obj.depthArray;
            
            Phi_0 = 0.75;
            Phi_inf = 0.05;
            B = 1600; % meters
            
            Rho_fluid = 1.024; % g/cc, seawater
            Rho_grain = 2.7;   % g/cc, smectite
            
            porosity = Phi_inf + (Phi_0 - Phi_inf)*exp(-depth./B);
            bulkDensity = porosity*Rho_fluid + (1 - porosity)*Rho_grain;
        end
        function [ pcgwInterp ] = CalcPcgw( obj , nonwettingSaturation )
            % returns pcgw in MPa
            
            pcgw1 = interp1( obj.MICP1.S_nw , obj.MICP1.Pc_gw , nonwettingSaturation );
            pcgw2 = interp1( obj.MICP2.S_nw , obj.MICP2.Pc_gw , nonwettingSaturation );
            
            pcgwInterp = mean( [pcgw1 pcgw2] , 2 );
        end
        function [ pchwInterp ] = CalcPchw( obj , nonwettingSaturation )
            % returns pchw in MPa

            pchw1 = interp1( obj.MICP1.S_nw , obj.MICP1.Pc_hw , nonwettingSaturation );
            pchw2 = interp1( obj.MICP2.S_nw , obj.MICP2.Pc_hw , nonwettingSaturation );
            
            pchwInterp = mean( [pchw1 pchw2] , 2 );
        end
        
        function [ radiusInMeters1 , radiusInMeters2 ] = CalcPoreRadiusFromPcgw( obj )
            % pcgw in MPa
            % radius in meters
            
            pcgw1 = obj.MICP1.Pc_gw;
            pcgw2 = obj.MICP2.Pc_gw;
            radiusInMicrons1 = 2 * 0.072 ./ pcgw1;
            radiusInMicrons2 = 2 * 0.072 ./ pcgw2;
            
            radiusInMeters1 = radiusInMicrons1 ./ 1e6;
            radiusInMeters2 = radiusInMicrons2 ./ 1e6;            
        end
        function [ slope1 , slope2 ] = CalcSlopeOfCumPSD( obj, stringType )
            % slope is in volume fraction per radius meter

            [radiusInMeters1, radiusInMeters2] = obj.CalcPoreRadiusFromPcgw();
            
            sNw1 = obj.MICP1.S_nw;
            sNw2 = obj.MICP2.S_nw;
            
            switch stringType
                case 'linear'
                    slope1 =     -( sNw1(2:end) - sNw1(1:end - 1) ) ...
                              ./ ( radiusInMeters1(2:end) - radiusInMeters1(1:end - 1) );
                case 'log'
                    slope1 =     -( sNw1(2:end) - sNw1(1:end - 1) ) ...
                              ./ ( log10(radiusInMeters1(2:end)) - log10(radiusInMeters1(1:end - 1)) );
            end
            switch stringType
                case 'linear'
                    slope2 =     -( sNw2(2:end) - sNw2(1:end - 1) ) ...
                              ./ ( radiusInMeters2(2:end) - radiusInMeters2(1:end - 1) );
                case 'log'
                    slope2 =     -( sNw2(2:end) - sNw2(1:end - 1) ) ...
                              ./ ( log10(radiusInMeters2(2:end)) - log10(radiusInMeters2(1:end - 1)) );
            end
        end
        
        %%% Plotting subclass functions
        function [ solFigure ] = PlotSol( obj , solFigure , exportTable )
            solFigure = obj.PlotSol@BCFormation( solFigure , exportTable );
            
            figure(solFigure)
            axis(obj.solAxis)
            title('Blake Ridge - Solubility Path')
        end
        function [ sat2PFigure ] = PlotSat2P( obj , sat2PFigure , exportTable , transitionZoneProperties , lineStyle )
            sat2PFigure = obj.PlotSat2P@BCFormation( sat2PFigure , exportTable , transitionZoneProperties , lineStyle );
            
            figure(sat2PFigure)
            axis(obj.satAxis)
            title('Blake Ridge - 2 Phase Case')
        end
        function [ sat3PFigure ] = PlotSat3P( obj , sat3PFigure , exportTable , lineStyle )
            sat3PFigure = PlotSat3P@BCFormation( obj , sat3PFigure , exportTable , lineStyle );
            
            figure(sat3PFigure)
            axis(obj.satAxis)
            title('Blake Ridge - 3 Phase Case')
        end
        
        function PlotMICP( obj )
            figure
            semilogy(1 - obj.MICP1.S_nw, obj.MICP1.Pc_gw, 'Linewidth', 2)
            hold on
            semilogy(1 - obj.MICP2.S_nw, obj.MICP2.Pc_gw, 'Linewidth', 2)
            xlabel('1 - S_n_w, or S_w')
            ylabel('Pc in MPa')
            title('Primary Drainage Capillary Pressure Curve')
        end
        function PlotCumPSD( obj )
            
            [radiusInMeters1, radiusInMeters2] = obj.CalcPoreRadiusFromPcgw();

            figure
            semilogx(radiusInMeters1, obj.MICP1.S_nw, 'Linewidth', 2)
            hold on
            semilogx(radiusInMeters2, obj.MICP2.S_nw, 'Linewidth', 2)
            xlabel('Pore radius in meters')
            ylabel('S_n_w')
            title('Cumulative Pore Size Distribution')
            axis([1e-9, 1e-4, 0, 1])
            set(gca, 'XDir', 'reverse')
            
        end
        function PlotPSD( obj , stringType )
            
            [radiusInMeters1, radiusInMeters2] = obj.CalcPoreRadiusFromPcgw();
            [slope1, slope2] = obj.CalcSlopeOfCumPSD( stringType );
            
            
            figure
            semilogx(radiusInMeters1(1:end - 1), slope1, 'Linewidth', 2)
            hold on
            semilogx(radiusInMeters2(1:end - 1), slope2, 'Linewidth', 2)
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
        
        
        function [ pcgwFigure ] = PlotPcgw( obj , pcgwFigure , exportTable , transitionZoneProperties , lineStylePc )
            pcgwFigure = PlotPcgw@BCFormation( obj , pcgwFigure , exportTable , transitionZoneProperties , lineStylePc );
            
            figure(pcgwFigure)

            axis(obj.pcgwAxis)
            title('Blake Ridge Gas Overpressure')
        end
        
        
        % unfinished plotting functions below
        function [ ratioFigure ] = PlotRatio( obj , ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio )
            [ ratioFigure ] = PlotRatio@BCFormation( obj , ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio );
            
            figure(ratioFigure)
                        
            title('Blake Ridge Overpressure Ratio')
            % legend('')
            axis(obj.ratioAxis)

            
            
        end
        function [ ratioFigure ] = PlotFractureRatio1( ~ , ratioFigure )
            figure(ratioFigure)
            
            plot( [1 1] , [300 600] , 'k--' , 'linewidth' , 2 )
            hold on
        end
    end
    methods (Static)
        function [ result ] = LoadLDT()
            
            % LDT.file_ID = fopen(uigetfile('*.dat','Select Lithodensity Log'));
            LDT.file_ID = fopen('164-995B_hldt.dat');
            
            LDT.Header1_HOLE = fgetl(LDT.file_ID);
            LDT.Header2_LEG = fgetl(LDT.file_ID);
            LDT.Header3_TOP = fgetl(LDT.file_ID);
            LDT.Header4_BOTTOM = fgetl(LDT.file_ID);
            
            % Add functionality to determine SizeCapture
            LDT.Header5_TOOLNAME = fgets(LDT.file_ID);
            LDT.Header6_TOOLUNITS = fgetl(LDT.file_ID);
            
            LDT.TOP_DEPTH = sscanf(LDT.Header3_TOP,'%*s %f');
            LDT.BOTTOM_DEPTH = sscanf(LDT.Header4_BOTTOM,'%*s %f');
            LDT.ScannedData = fscanf(LDT.file_ID,'%f',[4 Inf])';
            fclose(LDT.file_ID);
            
            LDT.Depth = LDT.ScannedData(:,1);
            LDT.RHOB = LDT.ScannedData(:,2);
            LDT.PEF = LDT.ScannedData(:,3);
            LDT.DRHO = LDT.ScannedData(:,4);
            
            LDT.Fixed_Depth_TOP = ceil(LDT.TOP_DEPTH);
            LDT.Fixed_Depth_BOTTOM = floor(LDT.BOTTOM_DEPTH);
            LDT.Depth_Interval = LDT.Fixed_Depth_BOTTOM - LDT.Fixed_Depth_TOP + 1;
            LDT.Fixed_Depth = linspace(LDT.Fixed_Depth_TOP,LDT.Fixed_Depth_BOTTOM,LDT.Depth_Interval)';
            LDT.Fixed_RHOB = interp1(LDT.Depth,LDT.RHOB,LDT.Fixed_Depth);
            LDT.Fixed_PEF = interp1(LDT.Depth,LDT.PEF,LDT.Fixed_Depth);
            LDT.Fixed_DRHO = interp1(LDT.Depth,LDT.DRHO,LDT.Fixed_Depth);
            LDT.Calc_DPHI_SS = (2.65-LDT.Fixed_RHOB)/(2.65-1);
            LDT.Calc_DPHI_LS = (2.71-LDT.Fixed_RHOB)/(2.71-1);
            
            result = LDT;
        end
        function [ result ] = LoadDIT()
            % DIT.file_ID = fopen(uigetfile('*.dat','Select Dual Induction Tool Log'));
            DIT.file_ID = fopen('164-995B_dit.dat');
            
            DIT.Header1_HOLE = fgetl(DIT.file_ID);
            DIT.Header2_LEG = fgetl(DIT.file_ID);
            DIT.Header3_TOP = fgetl(DIT.file_ID);
            DIT.Header4_BOTTOM = fgetl(DIT.file_ID);
            
            % Add functionality to determine SizeCapture
            DIT.Header5_TOOLNAME = fgets(DIT.file_ID);
            DIT.Header6_TOOLUNITS = fgetl(DIT.file_ID);
            
            DIT.TOP_DEPTH = sscanf(DIT.Header3_TOP,'%*s %f');
            DIT.BOTTOM_DEPTH = sscanf(DIT.Header4_BOTTOM,'%*s %f');
            DIT.ScannedData = fscanf(DIT.file_ID,'%f',[4 Inf])';
            fclose(DIT.file_ID);
            
            DIT.Depth = DIT.ScannedData(:,1);
            DIT.IDPH = DIT.ScannedData(:,2);
            DIT.IMPH = DIT.ScannedData(:,3);
            DIT.SFLU = DIT.ScannedData(:,4);
            
            DIT.Fixed_Depth_TOP = ceil(DIT.TOP_DEPTH);
            DIT.Fixed_Depth_BOTTOM = floor(DIT.BOTTOM_DEPTH);
            DIT.Depth_Interval = DIT.Fixed_Depth_BOTTOM - DIT.Fixed_Depth_TOP + 1;
            DIT.Fixed_Depth = linspace(DIT.Fixed_Depth_TOP,DIT.Fixed_Depth_BOTTOM,DIT.Depth_Interval)';
            DIT.Fixed_IDPH = interp1(DIT.Depth,DIT.IDPH,DIT.Fixed_Depth);
            DIT.Fixed_IMPH = interp1(DIT.Depth,DIT.IMPH,DIT.Fixed_Depth);
            DIT.Fixed_SFLU = interp1(DIT.Depth,DIT.SFLU,DIT.Fixed_Depth);
            
            result = DIT;
        end
        function [ result ] = LoadGR()
            % GR.file_ID = fopen(uigetfile('*.dat','Select Gamma Ray Log'));
            GR.file_ID = fopen('164-995B_ngt-dit.dat');

            GR.Header1_HOLE = fgetl(GR.file_ID);
            GR.Header2_LEG = fgetl(GR.file_ID);
            GR.Header3_TOP = fgetl(GR.file_ID);
            GR.Header4_BOTTOM = fgetl(GR.file_ID);

            % Add functionality to determine SizeCapture
            GR.Header5_TOOLNAME = fgets(GR.file_ID);
            GR.Header6_TOOLUNITS = fgetl(GR.file_ID);

            GR.TOP_DEPTH = sscanf(GR.Header3_TOP,'%*s %f');
            GR.BOTTOM_DEPTH = sscanf(GR.Header4_BOTTOM,'%*s %f');
            GR.ScannedData = fscanf(GR.file_ID,'%f',[11 Inf])';
            fclose(GR.file_ID);

            GR.Depth = GR.ScannedData(:,1);
            GR.SGR = GR.ScannedData(:,2);
            GR.CGR = GR.ScannedData(:,3);
            GR.THOR = GR.ScannedData(:,4);
            GR.URAN = GR.ScannedData(:,5);
            GR.POTA = GR.ScannedData(:,6);
            GR.W1NG = GR.ScannedData(:,7);
            GR.W2NG = GR.ScannedData(:,8);
            GR.W3NG = GR.ScannedData(:,9);
            GR.W4NG = GR.ScannedData(:,10);
            GR.W5NG = GR.ScannedData(:,11);

            GR.Fixed_Depth_TOP = ceil(GR.TOP_DEPTH);
            GR.Fixed_Depth_BOTTOM = floor(GR.BOTTOM_DEPTH);
            GR.Depth_Interval = GR.Fixed_Depth_BOTTOM - GR.Fixed_Depth_TOP + 1;
            GR.Fixed_Depth = linspace(GR.Fixed_Depth_TOP,GR.Fixed_Depth_BOTTOM,GR.Depth_Interval)';

            GR.Fixed_SGR = interp1(GR.Depth,GR.SGR,GR.Fixed_Depth);
            GR.Fixed_CGR = interp1(GR.Depth,GR.CGR,GR.Fixed_Depth);
            GR.Fixed_THOR = interp1(GR.Depth,GR.THOR,GR.Fixed_Depth);
            GR.Fixed_URAN = interp1(GR.Depth,GR.URAN,GR.Fixed_Depth);
            GR.Fixed_POTA = interp1(GR.Depth,GR.POTA,GR.Fixed_Depth);
            GR.Fixed_W1NG = interp1(GR.Depth,GR.W1NG,GR.Fixed_Depth);
            GR.Fixed_W2NG = interp1(GR.Depth,GR.W2NG,GR.Fixed_Depth);
            GR.Fixed_W3NG = interp1(GR.Depth,GR.W3NG,GR.Fixed_Depth);
            GR.Fixed_W4NG = interp1(GR.Depth,GR.W4NG,GR.Fixed_Depth);
            GR.Fixed_W5NG = interp1(GR.Depth,GR.W5NG,GR.Fixed_Depth);
            
            result = GR;
        end
        function [ result ] = LoadBRG()
            % BRG.file_ID = fopen(uigetfile('*.dat','Select Sonic Log'));
            BRG.file_ID = fopen('164-995B_sst-brg.dat');

            BRG.Header1_TITLE = fgetl(BRG.file_ID);
            BRG.Header2_HOLE = fgetl(BRG.file_ID);
            BRG.Header3_LEG = fgetl(BRG.file_ID);
            BRG.Header4_TOOL = fgetl(BRG.file_ID);
            BRG.Header5_OPERATOR = fgetl(BRG.file_ID);
            BRG.Header6_DATE = fgetl(BRG.file_ID);
            BRG.Header7_BLANK = fgetl(BRG.file_ID);
            BRG.Header8_NOTE = fgetl(BRG.file_ID);
            BRG.Header9_DATA = fgetl(BRG.file_ID);
            BRG.Header10_BLANK = fgetl(BRG.file_ID);
            BRG.Header11_UNITS = fgetl(BRG.file_ID);

            BRG.ScannedData = fscanf(BRG.file_ID,'%f',[3 Inf])';

            BRG.Depth = flip(BRG.ScannedData(:,1));
            BRG.VP = flip(BRG.ScannedData(:,2));
            BRG.VS = flip(BRG.ScannedData(:,3));

            BRG.TOP_DEPTH = BRG.Depth(1);
            BRG.BOTTOM_DEPTH = BRG.Depth(size(BRG.Depth,1));

            fclose(BRG.file_ID);

            BRG.Fixed_Depth_TOP = ceil(BRG.TOP_DEPTH);
            BRG.Fixed_Depth_BOTTOM = floor(BRG.BOTTOM_DEPTH);
            BRG.Depth_Interval = BRG.Fixed_Depth_BOTTOM - BRG.Fixed_Depth_TOP + 1;
            BRG.Fixed_Depth = linspace(BRG.Fixed_Depth_TOP,BRG.Fixed_Depth_BOTTOM,BRG.Depth_Interval)';

            BRG.Fixed_VP = interp1(BRG.Depth,BRG.VP,BRG.Fixed_Depth);
            BRG.Fixed_VS = interp1(BRG.Depth,BRG.VS,BRG.Fixed_Depth);
            
            result = BRG;
        end
        function [ CALI ] = LoadCALIHLDT()
            % DIT.file_ID = fopen(uigetfile('*.dat','Select Dual Induction Tool Log'));
            CALI.file_ID = fopen('164-995B_cali-hldt.dat');

            CALI.Header1_HOLE = fgetl(CALI.file_ID);
            CALI.Header2_LEG = fgetl(CALI.file_ID);
            CALI.Header3_TOP = fgetl(CALI.file_ID);
            CALI.Header4_BOTTOM = fgetl(CALI.file_ID);
            CALI.Header5_TOOLNAME = fgets(CALI.file_ID);
            CALI.Header6_TOOLUNITS = fgetl(CALI.file_ID);

            CALI.TOP_DEPTH = sscanf(CALI.Header3_TOP,'%*s %f');
            CALI.BOTTOM_DEPTH = sscanf(CALI.Header4_BOTTOM,'%*s %f');
            CALI.ScannedData = fscanf(CALI.file_ID,'%f',[2 Inf])';
            fclose(CALI.file_ID);

            CALI.Depth = CALI.ScannedData(:,1);
            CALI.Diameter = CALI.ScannedData(:,2);


            % CALI.Fixed_Depth_TOP = ceil(CALI.TOP_DEPTH);
            % CALI.Fixed_Depth_BOTTOM = floor(CALI.BOTTOM_DEPTH);
            % CALI.Depth_Interval = CALI.Fixed_Depth_BOTTOM - CALI.Fixed_Depth_TOP + 1;
            % CALI.Fixed_Depth = linspace(CALI.Fixed_Depth_TOP,CALI.Fixed_Depth_BOTTOM,CALI.Depth_Interval)';
            % CALI.Fixed_IDPH = interp1(CALI.Depth,CALI.IDPH,CALI.Fixed_Depth);
            % CALI.Fixed_IMPH = interp1(CALI.Depth,CALI.IMPH,CALI.Fixed_Depth);
            % CALI.Fixed_SFLU = interp1(CALI.Depth,CALI.SFLU,CALI.Fixed_Depth);
        end
        function [ result ] = LoadMICP1()
            MICP1.file_ID = fopen('MICP1_BR.dat');
            MICP1.ScannedData = fscanf(MICP1.file_ID,'%f',[3 Inf])';
            fclose(MICP1.file_ID);
            
            MICP1.S_nw = MICP1.ScannedData(:,1);
            MICP1.Pc_gw = MICP1.ScannedData(:,2);
            MICP1.Pc_hw = MICP1.ScannedData(:,3);
            
            result = MICP1;
        end
        function [ result ] = LoadMICP2()
            MICP2.file_ID = fopen('MICP2_BR.dat');
            MICP2.ScannedData = fscanf(MICP2.file_ID,'%f',[3 Inf])';
            fclose(MICP2.file_ID);
            
            MICP2.S_nw = MICP2.ScannedData(:,1);
            MICP2.Pc_gw = MICP2.ScannedData(:,2);
            MICP2.Pc_hw = MICP2.ScannedData(:,3);
            
            result = MICP2;
        end
    end
    % UNUSED CLASS METHODS
    %{
        function LoadSolubility( obj )
            load('Blake Ridge Data\Solubility Plots\BR_bulk_C_L_G.mat')
            obj.Bulk.Solubility = BR_bulk_C_L_G(:,1); % mol CH4/kg H2O
            obj.Bulk.Depth = BR_bulk_C_L_G(:,2) - obj.seafloorDepth; % mbsf
                        
        end        
        function PreprocessData( obj )
            
            obj.Data.depth_top = min([obj.LDT.Fixed_Depth_TOP   obj.Depth.saturationTop]);
            obj.Data.depth_bottom = max([obj.LDT.Fixed_Depth_BOTTOM obj.Depth.saturationBottom]);
            
            obj.Saturation.Above_Data = repmat( obj.Saturation.HydrateGrid(1,:) , obj.Depth.Graph_Top - obj.Depth.saturationTop , 1 );
            obj.Saturation.Below_Data = repmat( obj.Saturation.GasGrid(end,:) , obj.Depth.saturationBottom - obj.Depth.Graph_Bottom , 1 );
            
            obj.Saturation.Hydrate_Data = vertcat(obj.Saturation.Above_Data,obj.Saturation.HydrateGrid);
            obj.Saturation.Gas_Data = vertcat(obj.Saturation.GasGrid,obj.Saturation.Below_Data);
            

            obj.Data.depth_interval = obj.Data.depth_bottom - obj.Data.depth_top + 1;
            obj.Data.log = zeros( obj.Data.depth_interval , 50 );
            obj.Data.log(:,1) = linspace( obj.Data.depth_top , obj.Data.depth_bottom , obj.Data.depth_interval );

%             obj.DataTable.depth = linspace( obj.Data.depth_top , obj.Data.depth_bottom , obj.Data.depth_interval )';
            obj.DataTable.depth = (obj.Data.depth_top : 1 : obj.Data.depth_bottom)';
        end
            
    %}
    % UNUSED STATIC METHODS
    %{
        function [ result ] = LoadPhaseBehaviorBlakeRidge()
            % Loads saturation data from hydrate and gas grabit array .mat files and
            %   saves into struct Saturation with fields 'Hydrate' and 'Gas'
            Saturation = struct('Hydrate',cell2mat(struct2cell(Load('Blake Ridge Data\Phase Behavior\HydrateBehavior.mat'))),'Gas',cell2mat(struct2cell(Load('Blake Ridge Data\Phase Behavior\GasBehavior.mat'))));
            
            % Upper and lower depths of grabit depth of investigation
            Depth_grid = linspace(3200,3300,101);
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
            
            % Hard coded fixes
            Saturation.GasGrid(49,40) = 0;
            
            Saturation.GasGrid(66:end,5) = Saturation.GasGrid(end,5);
            Saturation.GasGrid(66:end,6) = Saturation.GasGrid(end,6);
            Saturation.GasGrid(65:end,7) = Saturation.GasGrid(end,7);
            Saturation.GasGrid(65:end,8) = Saturation.GasGrid(end,8);
            Saturation.GasGrid(65:end,9) = Saturation.GasGrid(end,9);
            Saturation.GasGrid(67:end,10) = Saturation.GasGrid(end,10);

            Saturation.Three_Phase_Top = zeros(1,size(Saturation.HydrateGrid,2));
            Saturation.Three_Phase_Bottom = zeros(1,size(Saturation.HydrateGrid,2));
            for index=1:size(Saturation.HydrateGrid,2)
                Saturation.Three_Phase_Top(index) = sum(~any(Saturation.GasGrid(:,index),2)) + 1;
                Saturation.Three_Phase_Bottom(index) = sum(any(Saturation.HydrateGrid(:,index),2));
            end
            
            result = Saturation;
            
            
        end
    %}
end































