classdef DCSeismicAnalysisBR < DCBlakeRidge
    properties
        graphTop
        graphBottom
        saturationBottom
        
        quantityArray
        
        SaturationLF
        Dickens
    end
    properties (Constant)
        liuTop = 3200; % mbsl
        liuBottom = 3300; % mbsl
        saturationTop = 200; % mbsf
        
        sandK = 38*10^9; % Pa
        clayK = 21.2*10^9; % Pa
        waterK = 2.688*10^9; % Pa
        
    end
    methods
        %%% Constructor
        function [ obj ] = DCSeismicAnalysisBR()
            obj@DCBlakeRidge()
            
            obj.graphTop = obj.liuTop - obj.seafloorDepth; % mbsf
            obj.graphBottom = obj.liuBottom - obj.seafloorDepth; % mbsf
            obj.saturationBottom = obj.maxDepth; % mbsf

            obj.depthArray = (obj.minDepth : 1 : obj.maxDepth)';
            
            obj.quantityArray = (1 : 1 : 40)';
            
            obj.SaturationLF = obj.LoadPhaseBehaviorBlakeRidge();
            
%             obj.Dickens = obj.LoadDickensBlakeRidge();
            
            
        end
        
        %%% Main function
        function RunSeismicAnalysisRoutine( obj )
            %%% Make main table
            data = table();
            data.Depth = obj.depthArray;
            nDepth = numel(obj.depthArray);
            
            data.Pressure = obj.CalcPressure(data.Depth);
            
            %%% Load log data into table
            data.Resistivity = nan(nDepth, 1);
            data.Resistivity(obj.DIT.Fixed_Depth_TOP : obj.DIT.Fixed_Depth_BOTTOM) = obj.DIT.Fixed_IDPH;
            
            data.GammaRay = nan(nDepth, 1);
            data.GammaRay(obj.GR.Fixed_Depth_TOP : obj.GR.Fixed_Depth_BOTTOM) = obj.GR.Fixed_SGR;
            
            % VP in m^3/s
            data.VP = nan(nDepth, 1);
            data.VP(obj.BRG.Fixed_Depth_TOP : obj.BRG.Fixed_Depth_BOTTOM) = obj.BRG.Fixed_VP;
            
            % VS in m^3/s
            data.VS = nan(nDepth, 1);
            data.VS(obj.BRG.Fixed_Depth_TOP : obj.BRG.Fixed_Depth_BOTTOM) = obj.BRG.Fixed_VS;
            
            % Density in kg/m^3
            data.BulkDensity = nan(nDepth, 1);
            data.BulkDensity(obj.LDT.Fixed_Depth_TOP : obj.LDT.Fixed_Depth_BOTTOM) = obj.LDT.Fixed_RHOB .* 1000;

            %%% Methane quantity for loop
            nQuantity = numel(obj.quantityArray);
            for iQuantity = 1:nQuantity
                quantity = obj.quantityArray(iQuantity);
                

                
                data.Sh = obj.SaturationLF.Hydrate(:, iQuantity);
                data.Sg = obj.SaturationLF.Gas(:, iQuantity);
                data.Sw = 1 - data.Sh - data.Sg;
                
                
                % Gets indices of depths at BSR and BGHSZ
                Wave.BSR = find(data.Depth == obj.graphTop + obj.SaturationLF.Top3P(iQuantity));
                Wave.BGHSZ = find(data.Depth == obj.graphTop + obj.SaturationLF.Bottom3P(iQuantity));
                Wave.rickerFrequency = 30;
                
                data.Porosity = obj.CalcPorosity(data.Resistivity, data.Sw);
                
                % Gas bulk modulus in Pa
                data.GasK = data.Pressure;
                
                % Bulk modulus from VP and VS in Pa
                data.BulkK = data.BulkDensity .* (data.VP .^ 2 - (4/3) .* (data.VS .^ 2));
                
                
                
                hydrateK = 6.414*10^9; % Pa
                
                % Fluid modulus using isostress average (volume-weighted harmonic mean) 
                data.FluidK = 1 ./ ...
                                (data.Sw ./ obj.waterK ...
                               + data.Sg ./ data.GasK ...
                               + data.Sh ./ hydrateK);
                
                data.GrainK = obj.CalcGrainK(data.GammaRay);
                
                data.FrameK = obj.CalcFrameK(data.Porosity, data.GrainK);
                
                data.ShearG = data.BulkDensity .* data.VS .^ 2;
                
                %%% Fluid substitution
                data.BulkKFS = obj.CalcBulkKFluidSubstituted(data.Porosity, data.FluidK, data.GrainK, data.FrameK);
                
                data.BulkDensityFS = obj.CalcBulkDensityFluidSubstituted(data.GammaRay, data.Porosity, data.Sw, data.Sg, data.Sh);
                
                data.VSFS = (data.ShearG ./ data.BulkDensityFS) .^ (1/2);
                
                data.VPFS = ( data.BulkKFS ./ data.BulkDensityFS ...
                                + (4/3) .* data.VSFS .^ 2 ...
                            ) .^ (1/2);
                
                Wave = obj.FindIntervalOfOneSeismicWavelength(data.VPFS, Wave);
                
                
            end
            
            
            
        end
        
        
        %%% Rock physics calculations
        function [ Porosity ] = CalcPorosity( obj , Rt , Sw )
            % Calculates Rw based on Ro, To, and T at the specified depth

            depth = obj.depthArray ./ 1000;         % convert to km
            temperatureGradient995 = 38.5;   % C deg/km (36.9 for well 997)
            seafloorTemperature995 = 3;  % C deg

            referenceRo = .24;         % ohm-m (.223)
            referenceTemp = 18;        % C deg, given in input geophysics file

            a = 0.9;
            m = 2.7;
            n = 1.9386;

            temperature = seafloorTemperature995 + depth .* temperatureGradient995;
            Rw = referenceRo .* (referenceTemp + 21.5) ./ (temperature + 21.5);
            Porosity = (a .* Rw ./ Rt ...
                            ./ (Sw .^ n)) ...
                            .^ (1 ./ m);
            
            %%% OLD CODE
            %{
            Assumed_PPM = 32;

            % best fit function
            Ro = .8495 + 2.986e-4.*Depth*1000;

            Saturation_Hydrate = 1 - (Ro./Resistivity_t).^(1/Coeff_n);
            Saturation_Hydrate(Saturation_Hydrate<0)=0;

            Old_Porosity = (Coeff_a.*Rw./Resistivity_t).^(1/Coeff_m);
            Old_Porosity_Corrected = Old_Porosity./(1 - Saturation_Hydrate);
            %}
        end
        function [ grainK ] = CalcGrainK( obj , gamma )  
            % Assumed normalizing terms
            shaleFreeSGR = 12;
            shalePureSGR = 110;
            
            distanceSGR = shalePureSGR - shaleFreeSGR;
            gammaNormalized = (gamma - shaleFreeSGR) ./ distanceSGR;
            gammaNormalized(gammaNormalized < 0) = 0;
            gammaNormalized(gammaNormalized > 1) = 1;
            
            voigtK = gammaNormalized .* obj.clayK ...
                    + (1 - gammaNormalized) .* obj.sandK;
            reussK = (obj.sandK .* obj.clayK) ...
                        ./ (gammaNormalized .* obj.sandK ...
                            + (1 - gammaNormalized) .* obj.clayK);
            grainK = (1/2) .* (voigtK + reussK);
        end
        function [ frameK ] = CalcFrameK( ~ , porosity , grainK )
            
            frameK = grainK ...
                    .* 10 .^ (3.02 - 7.372 .* porosity);
            return
            
            frameK = grainK .* (bulkK .* (porosity .* (grainK - fluidK) + fluidK) - grainK .* fluidK) ...
                ./ (porosity .* grainK .* (grainK - fluidK) + fluidK .* (bulkK - grainK));
        end
        function [ bulkKFS ] = CalcBulkKFluidSubstituted( ~ , porosity , fluidK , grainK , frameK )
            
            Q = (fluidK .* (grainK - frameK)) ...
                ./ (porosity .* (grainK - fluidK));
            
            bulkKFS = grainK .* ((frameK + Q) ...
                                ./ (grainK + Q));   
        end
        function [ bulkDensityFS ] = CalcBulkDensityFluidSubstituted( obj , gamma , porosity , Sw , Sg , Sh )
            waterDensity = obj.waterDensity / 1000; % g/cm^3
            gasDensity = 0.3; % g/cm^3
            hydrateDensity = (0.924 + 0.933) / 2; % g/cm^3
            clayDensity = 2.6; % g/cm^3
            sandDensity = 2.7; % g/cm^3
            
            shaleFreeSGR = 12;
            shalePureSGR = 110;
            
            distanceSGR = shalePureSGR - shaleFreeSGR;
            gammaNormalized = (gamma - shaleFreeSGR)./distanceSGR;
            gammaNormalized(gammaNormalized < 0) = 0;
            gammaNormalized(gammaNormalized > 1) = 1;
            
            porosity(porosity < 0) = 0;
            porosity(porosity > 1) = 1;
            
            fluidDensity = Sw .* waterDensity ...
                         + Sg .* gasDensity ...
                         + Sh .* hydrateDensity;
             
            matrixDensity = gammaNormalized .* clayDensity ...
                            + (1 - gammaNormalized) .* sandDensity;
            
            bulkDensityFS = porosity .* fluidDensity ...
                            + (1 - porosity) .* matrixDensity;
            
            bulkDensityFS = bulkDensityFS .* 1000; % convert g/cm^3 to kg/m^3
        end
        function [ Wave ] = FindIntervalOfOneSeismicWavelength( ~ , VPFS , Wave )
            % This finds the index above the BSR and the index below the
            % BGHSZ to have an interval over which we average the fluid
            % substituted parameters
            % 
            % Calculates the average VP starting with depth one meter above
            % the BSR and finds the wavelength using lambda = VP/freqency
            %
            % Iterates by increasing this search interval until the
            % calculated wavelength lambda = the thickness of the search
            % interval
            
            if isempty(Wave.BSR) || isempty(Wave.BGHSZ)
                Wave.topWavelengthIndex = [];
                Wave.bottomWavelengthIndex = [];
                return
            end
            
            
            iterateTop = true;
            iterateBottom = true;
            iterationThickness = 0;
            
            while iterateTop || iterateBottom
                iterationThickness = iterationThickness + 1;
                
                topVelocity = mean(VPFS(Wave.BSR - iterationThickness : Wave.BSR - 1));
                topWavelength = topVelocity / Wave.rickerFrequency;
                
                bottomVelocity = mean(VPFS(Wave.BGHSZ + 1 : Wave.BGHSZ + iterationThickness));
                bottomWavelength = bottomVelocity / Wave.rickerFrequency;
                
                if iterateTop && (topWavelength - iterationThickness) < 0.5
                    Wave.topWavelengthIndex = Wave.BSR - iterationThickness;
                    iterateTop = false;
                end
                
                if iterateBottom && (bottomWavelength - iterationThickness) < 0.5
                    Wave.bottomWavelengthIndex = Wave.BGHSZ + iterationThickness;
                    iterateBottom = false;
                end
            end
        end
        
        %%% Loading methods
        function [ SaturationLF ] = LoadPhaseBehaviorBlakeRidge( obj )
            % Loads saturation data from hydrate and gas grabit array .mat files and
            % saves into struct Saturation with fields 'Hydrate' and 'Gas'
            s = struct('Hydrate', cell2mat(struct2cell(load('Blake Ridge Data\Phase Behavior\HydrateBehavior.mat'))), ...
                       'Gas', cell2mat(struct2cell(load('Blake Ridge Data\Phase Behavior\GasBehavior.mat'))));

            % Upper and lower depths of grabit depth of investigation
            depthGrid = linspace(obj.liuTop, obj.liuBottom, 101);
            methaneQuantityGrid = linspace(0, 40, 41);
            
            [X, Y] = meshgrid(methaneQuantityGrid, depthGrid);
            
            % Hydrate saturations
            HydrateGrid = griddata(s.Hydrate(:,1), s.Hydrate(:,2), s.Hydrate(:,3), X, Y, 'cubic');
            HydrateGrid(HydrateGrid < 0) = 0;
            HydrateGrid(isnan(HydrateGrid)) = 0;
            HydrateGrid(HydrateGrid < .0001 & HydrateGrid ~= 0) = 0;
            
            % Gas saturations
            GasGrid = griddata(s.Gas(:,1), s.Gas(:,2), s.Gas(:,3), X, Y, 'cubic');
            GasGrid(GasGrid < 0) = 0;
            GasGrid(isnan(GasGrid)) = 0;
            GasGrid(GasGrid < .0001 & GasGrid ~= 0) = 0;
            
            % Hard coded fixes
            GasGrid(49, 40) = 0;
            GasGrid(66:end, 5) = GasGrid(end, 5);
            GasGrid(66:end, 6) = GasGrid(end, 6);
            GasGrid(65:end, 7) = GasGrid(end, 7);
            GasGrid(65:end, 8) = GasGrid(end, 8);
            GasGrid(65:end, 9) = GasGrid(end, 9);
            GasGrid(67:end, 10) = GasGrid(end, 10);
            
            n = numel(obj.quantityArray);
            
            ThreePhaseTop = nan(1, n);
            ThreePhaseBottom = nan(1, n);
            
            HydrateFull = zeros(numel(obj.depthArray), n);
            GasFull = zeros(numel(obj.depthArray), n);
            
            for i = 1:n
                iQuantity = obj.quantityArray(i);
                
                hasGasLogical = GasGrid(:, iQuantity + 1) > 0;
                if any(hasGasLogical)
                    ThreePhaseTop(i) = find(hasGasLogical, 1, 'first');
                end
                
                hasHydrateLogical = HydrateGrid(:, iQuantity + 1) > 0;
                if any(hasHydrateLogical)
                    ThreePhaseBottom(i) = find(hasHydrateLogical, 1, 'last');
                end
                
                HydrateFull(obj.saturationTop:obj.graphTop - 1, i) = HydrateGrid(1, i);
                HydrateFull(obj.graphTop:obj.graphBottom, i) = HydrateGrid(:, i);
                
                GasFull(obj.graphTop:obj.graphBottom, i) = GasGrid(:, i);
                GasFull(obj.graphBottom + 1:obj.saturationBottom, i) = GasGrid(end, i);        
            end
            
            
            
            SaturationLF.Top3P = ThreePhaseTop;
            SaturationLF.Bottom3P = ThreePhaseBottom;
            
            SaturationLF.Hydrate = HydrateFull;
            SaturationLF.Gas = GasFull;
            
            
            
        end
        function [ Dickens ] = LoadDickensBlakeRidge( obj )
            % Returns CH4 quantity vs depth using Dickens et al. 1997
            % Input array is in mbsf
            % Output array is in rounded 

            % From grabit
            Dickens = load('Blake Ridge Data\Phase Behavior\MethaneQuantity.mat');
            % Column 1 is CH4 quantity
            % Colume 2 is mbsf
            Dickens.MethaneQuantity(:, 1) = Dickens.MethaneQuantity(:, 1) .* 16.04; % conversion from mol CH4 to g CH4 (per dm^3)

            % X is CH4 quantity
            % Y is mbsf
            [X, Y] = meshgrid(linspace(0, 40, 41), obj.depthArray);
            
            hydrateSat = [zeros(numel(obj.depthArray), 1), obj.SaturationLF.Hydrate];
            gasSat = [zeros(numel(obj.depthArray), 1), obj.SaturationLF.Gas];
            
            interpSh = interp2(X, Y, hydrateSat, Dickens.MethaneQuantity(:, 1), Dickens.MethaneQuantity(:, 2));
            interpSg = interp2(X, Y, gasSat, Dickens.MethaneQuantity(:, 1), Dickens.MethaneQuantity(:, 2));

            Dickens.Hydrate = interp1( Dickens.MethaneQuantity(:,2) , interpSh , obj.depthArray );
            Dickens.Gas = interp1( Dickens.MethaneQuantity(:,2) , interpSg , obj.depthArray );

            Dickens.Hydrate(isnan(Dickens.Hydrate)) = 0;
            Dickens.Gas(isnan(Dickens.Gas)) = 0;

            % finds average methane quantity over upper and lower bound of Dickens
            % upper_bound = 130;
            % lower_bound = 590;
            % interp_methane_quantity = interp1( Dickens.MethaneQuantity(:,2) , Dickens.MethaneQuantity(:,1) , depth_array( depth_array <= lower_bound & depth_array >= upper_bound ) );
            % interp_methane_quantity(isnan(interp_methane_quantity)) = [];
            % average_methane_quantity = mean(interp_methane_quantity);
            
            % Methane quantity vs depth
            figure

            axis_1 = subplot(2,1,1);
            h = plot( Dickens.MethaneQuantity(:,1) , Dickens.MethaneQuantity(:,2) , 'k' );
            set(gca,'YDir','Reverse')
            axis([0 40 100 600])
            % grid on
            xlabel('Methane quantity (g/dm^3 of pore volume)')
            ylabel('Depth (mbsf)')
            title('a) Experimentally Measured Methane Quantity')
            set(h,'LineWidth',2.5);

            % Fluid saturations vs depth
            % figure
            % scatter( Sh_interpolated , Dickens.MethaneQuantity(:,2) )
            % hold on
            % scatter( Sg_interpolated , Dickens.MethaneQuantity(:,2) )
            % hold on
            axis_2 = subplot(2,1,2);
            h = plot( Dickens.Hydrate , obj.depthArray );
            set(h,'LineWidth',2.5);
            set(h,'Color',[0 .5 0]);
            hold on
            h = plot( Dickens.Gas , obj.depthArray , 'r--' );
            set(h,'LineWidth',2.5);
            set(gca,'YDir','Reverse')
            axis([0 0.2 100 600])
            % axis([0 0.2 420 520])
            xlabel('Fluid saturations')
            ylabel('Depth (mbsf)')
            title('b) Variable Methane Quantity Scenario')
            % grid on
            legend('Hydrate saturation' , 'Gas saturation' )


            axis_1.Position = [.15 .59 .78 .37];
            axis_2.Position = [.15 .08 .78 .37];


            figure_1 = gcf;
            figure_1.Position(3) = 320;
            set(findall(figure_1,'-property','FontSize'),'FontSize',8)
            set(findall(figure_1,'-property','FontName'),'FontName','Arial')
        end        
        
        
        %%% Unused methods
        
        
        
    end
    methods (Static)
        
    end
end























