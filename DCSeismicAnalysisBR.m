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
        temperatureGradient995 = 38.5;   % C deg/km (36.9 for well 997)
        seafloorTemperature995 = 3;  % C deg
        
        liuTop = 3200; % mbsl
        liuBottom = 3300; % mbsl
        saturationTop = 200; % mbsf
        
        sandK = 38*10^9; % Pa
        clayK = 21.2*10^9; % Pa
        waterK = 2.688*10^9; % Pa
        
        selectedQuantities = [6 15 23 32 40];
        
        aHelgerudK = -1.09e-2; % GPa/C deg
        bHelgerudK = 3.8e-3; % GPa/MPa
        cHelgerudK = 8.39; % GPa
        
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
        
        %%% Main methods
        function [ Wave , data , WaveBase , dataBase ] = RunSeismicAnalysisRoutine( obj , caseString , baseFlag )
            %%% Instantiate optional base (Sw = 1) calculation variables
            WaveBase = struct();
            dataBase = table();
            
            
            %%% Instantiate Wave struct
            Wave.rickerFrequency = 30;

            
            
            %%% Instantiate main table
            data = table();
            data.Depth = obj.depthArray;
            nDepth = numel(obj.depthArray);
            
            
            
            data.Pressure = obj.CalcPressure(data.Depth);
            data.Temperature = obj.CalcTemperature995();
            
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
            
            % Density in kg/m^3 (converted from g/cm^3 in the log)
            data.BulkDensity = nan(nDepth, 1);
            data.BulkDensity(obj.LDT.Fixed_Depth_TOP : obj.LDT.Fixed_Depth_BOTTOM) = obj.LDT.Fixed_RHOB .* 1000;

            
            data.HydrateK = obj.CalcHydrateK(data.Pressure, data.Temperature);
            
            
%             newPorosityFlag = true;
            newPorosityFlag = false;
            if newPorosityFlag
                data.Porosity = obj.CalcPorosity(data.Resistivity, 1, data.Temperature);
            end
            
            %%% temporary code to get porosity for Sw = 1 for background
            %%% properties figure
            dataBase = data;
            dataBase.Porosity = obj.CalcPorosity(data.Resistivity, 1, data.Temperature);
            %%% end of temp code
            
            
            %%% Adjust for parameter sensitivity
            switch caseString
                case 'ParameterSensitivity'
                    upperBoundIndex = 415;
                    lowerBoundIndex = 630;
                    data.Resistivity(:) = mean(data.Resistivity(upperBoundIndex : lowerBoundIndex));
                    data.GammaRay(:) = mean(data.GammaRay(upperBoundIndex : lowerBoundIndex));
                    data.VP(:) = mean(data.VP(upperBoundIndex : lowerBoundIndex));
                    data.VS(:) = mean(data.VS(upperBoundIndex : lowerBoundIndex));
                    data.BulkDensity(:) = mean(data.BulkDensity(upperBoundIndex : lowerBoundIndex));
                    
                    data.HydrateK(:) = mean(data.HydrateK(upperBoundIndex : lowerBoundIndex));
            end
            
            
            
            %%% Methane quantity for loop
            nQuantity = numel(obj.quantityArray);
            
            % Instantiate storage variables
            Wave.seismogram = cell(nQuantity, 1);
            Wave.thickness3P = cell(nQuantity, 1);
            Wave.peak = cell(nQuantity, 2);
            Wave.time = cell(nQuantity, 1);
            Wave.depth = data.Depth;
            
            Wave.bulkKFS = cell(nQuantity, 1);
            Wave.bulkDensityFS = cell(nQuantity, 1);
            Wave.VSFS = cell(nQuantity, 1);
            Wave.VPFS = cell(nQuantity, 1);
            
            
            %%% Run base case for parameter sensitivity
            if baseFlag
                [ ~ , ~ , ~ , WaveBase ] ...
                        = obj.RunRockPhysicsRoutine(data, Wave, 0, caseString, newPorosityFlag);
            end
            
            %%% For loop through all methane quantities
            for iQuantity = 1:nQuantity
                quantity = obj.quantityArray(iQuantity);
                
                [ seismogram , timeSeries , thickness , parameterSensitivity ] ...
                    = obj.RunRockPhysicsRoutine(data, Wave, iQuantity, caseString, newPorosityFlag);
                
                
                Wave.seismogram{iQuantity} = seismogram;
                Wave.time{iQuantity} = timeSeries;
                Wave.thickness3P{iQuantity} = thickness;
                
                [ peakValue , peakIndex ] = findpeaks(seismogram);
                
                switch numel(peakValue)
                    case 1
                        Wave.peak{iQuantity, 1} = peakValue(1);
                    case 2
                        Wave.peak{iQuantity, 1} = peakValue(1);
                        Wave.peak{iQuantity, 2} = peakValue(2);
                end
                
                Wave.bulkKFS{iQuantity} = parameterSensitivity.bulkKFS;
                Wave.bulkDensityFS{iQuantity} = parameterSensitivity.bulkDensityFS;
                Wave.VSFS{iQuantity} = parameterSensitivity.VSFS;
                Wave.VPFS{iQuantity} = parameterSensitivity.VPFS;
                
                
            end
        end
        function [ seismogram , timeSeries , thickness , parameterSensitivity ] = ...
                RunRockPhysicsRoutine( obj , data , Wave , iQuantity , caseString , newPorosityFlag )
            
            seismogram = [];
            timeSeries = [];
            thickness = [];
            parameterSensitivity = [];
            
            
            Wave.BSR = [];
            Wave.BGHSZ = [];
            
            if iQuantity == 0
                data.Sh = zeros(height(data), 1);
                data.Sg = zeros(height(data), 1);
                data.Sw = ones(height(data), 1);
            else
                data.Sh = obj.SaturationLF.Hydrate(:, iQuantity);
                data.Sg = obj.SaturationLF.Gas(:, iQuantity);
                data.Sw = 1 - data.Sh - data.Sg;
                
                % Gets indices of depths at BSR and BGHSZ
                Wave.BSR = find(data.Depth == obj.graphTop + obj.SaturationLF.Top3P(iQuantity));
                Wave.BGHSZ = find(data.Depth == obj.graphTop + obj.SaturationLF.Bottom3P(iQuantity));
            end
            
            
            
            
            if ~newPorosityFlag
                data.Porosity = obj.CalcPorosity(data.Resistivity, data.Sw, data.Temperature);
            end
            
            % Gas bulk modulus in Pa
            data.GasK = obj.CalcGasK(data.Pressure);
            
            switch caseString
                case 'ParameterSensitivity'
                    upperBoundIndex = 415;
                    lowerBoundIndex = 630;
                    data.Porosity(:) = mean(data.Porosity(upperBoundIndex : lowerBoundIndex));
                    data.GasK(:) = mean(data.GasK(upperBoundIndex : lowerBoundIndex));
            end
            
            
            % Bulk modulus from VP and VS in Pa
            data.BulkK = obj.CalcBulkK(data.BulkDensity, data.VP, data.VS);
            hydrateK = 6.414*10^9; % Pa
            % Fluid modulus using isostress average (volume-weighted harmonic mean) 
            data.FluidK = obj.CalcFluidK(data.Sw, obj.waterK, data.Sg, data.GasK, data.Sh, data.HydrateK);
            data.GrainK = obj.CalcGrainK(data.GammaRay);
            data.FrameK = obj.CalcFrameK(data.Porosity, data.GrainK);
            data.ShearG = obj.CalcShearG(data.BulkDensity, data.VS);

            %%% Fluid substitution
            data.BulkKFS = obj.CalcBulkKFluidSubstituted(data.Porosity, data.FluidK, data.GrainK, data.FrameK);
            data.BulkDensityFS = obj.CalcBulkDensityFluidSubstituted(data.GammaRay, data.Porosity, data.Sw, data.Sg, data.Sh);
            data.VSFS = obj.CalcVSFS(data.ShearG, data.BulkDensityFS);
            data.VPFS = obj.CalcVPFS(data.BulkKFS, data.BulkDensityFS, data.VSFS);
            
            Wave = obj.FindIntervalOfOneSeismicWavelength(data.VPFS, Wave);
            originalResolutionFlag = true;
            [ data , Wave ] = obj.CalcAverageProperties(data, Wave, originalResolutionFlag);
            data.Impedance = obj.CalcImpedance(data.BulkDensityFSAvg, data.VPFSAvg);
            data.ReflectionCoefficient = obj.CalcReflectionCoefficient(data.Impedance);
            data.ConvolutionDt = obj.CalcConvolutionDt(data.VPFS, data.VPFSAvg, originalResolutionFlag);
            data.TimeSeries = obj.CalcTimeSeries(data.ConvolutionDt);
            
            [ ~ , ~ , seismogram ] = obj.CalcSeismogram( Wave.rickerFrequency , data.TimeSeries , data.ReflectionCoefficient );
            timeSeries = data.TimeSeries;
            thickness = data.Depth(Wave.BGHSZ) - data.Depth(Wave.BSR);
            
            switch caseString
                case 'ParameterSensitivity'
                    parameterSensitivity.bulkKFS = data.BulkKFS;
                    parameterSensitivity.bulkDensityFS = data.BulkDensityFS;
                    parameterSensitivity.VSFS = data.VSFS;
                    parameterSensitivity.VPFS = data.VPFS;
            end
            
            
        end
        
        %%% Rock physics calculation methods
        function [ temperature ] = CalcTemperature995( obj )
            depth = obj.depthArray ./ 1000; % convert to km
            temperature = obj.seafloorTemperature995 + depth .* obj.temperatureGradient995; % in C deg
        end
        function [ Porosity ] = CalcPorosity( obj , Rt , Sw , temperature )
            % Calculates Rw based on Ro, To, and T at the specified depth
%             depth = obj.depthArray ./ 1000;         % convert to km
            
            referenceRo = .24;         % ohm-m (.223)
            referenceTemp = 18;        % C deg, given in input geophysics file
            
            a = 0.9;
            m = 2.7;
            n = 1.9386;
            
%             temperature = obj.seafloorTemperature995 + depth .* obj.temperatureGradient995;
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
        function [ gasK ] = CalcGasK( ~ , hydrostaticPressure )
            gasK = hydrostaticPressure;
        end
        function [ hydrateK ] = CalcHydrateK( obj , hydrostaticPressure , temperature )
            a = obj.aHelgerudK;
            b = obj.bHelgerudK;
            c = obj.cHelgerudK;
            
            hydrateK = a .* temperature ...
                     + b .* hydrostaticPressure ./ 1e6 ...
                     + c; % in GPa
            hydrateK = hydrateK .* 1e9; % in Pa
        end
        function [ bulkK ] = CalcBulkK( ~ , bulkDensity , VP , VS )
            bulkK = bulkDensity .* (VP .^ 2 - (4/3) .* (VS .^ 2));
        end
        function [ fluidK ] = CalcFluidK( ~ , Sw , waterK , Sg , gasK , Sh , hydrateK )
            fluidK = 1 ./ ...
                    (Sw ./ waterK ...
                   + Sg ./ gasK ...
                   + Sh ./ hydrateK);
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
        function [ shearG ] = CalcShearG( ~ , bulkDensity , VS )
            shearG = bulkDensity .* VS .^ 2;
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
%             hydrateDensity = 0.9; % g/cm^3
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
        function [ VSFS ] = CalcVSFS( ~ , shearG , bulkDensityFS )
            VSFS = (shearG ./ bulkDensityFS) .^ (1/2);
        end
        function [ VPFS ] = CalcVPFS( ~ , bulkKFS , bulkDensityFS , VSFS )
            VPFS = ( bulkKFS ./ bulkDensityFS ...
                        + (4/3) .* VSFS .^ 2 ...
                   ) .^ (1/2);
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
        function [ data , Wave ] = CalcAverageProperties( obj , data , Wave , originalResolutionFlag )
            
            data.VPFSAvg = data.VPFS;
            data.BulkDensityFSAvg = data.BulkDensityFS;
            data.VSFSAvg = data.VSFS;
            
            if isempty(Wave.topWavelengthIndex) || isempty(Wave.bottomWavelengthIndex)
                return
            end
            
            Wave.topLimit = obj.saturationTop;
            
            % Averaged VP
            AVO.vp1 = mean(data.VPFSAvg(Wave.topWavelengthIndex : Wave.BSR - 1));
            AVO.vp2 = mean(data.VPFSAvg(Wave.BSR : Wave.BGHSZ));
            AVO.vp3 = mean(data.VPFSAvg(Wave.BGHSZ + 1 : Wave.bottomWavelengthIndex));
            
            data.VPFSAvg(Wave.topLimit : Wave.BSR - 1) = AVO.vp1;
            if ~originalResolutionFlag
                data.VPFSAvg(Wave.BSR : Wave.BGHSZ) = AVO.vp2;
            end
            data.VPFSAvg(Wave.BGHSZ + 1 : end) = AVO.vp3;
            
            
            % Averaged Density
            AVO.density1 = mean(data.BulkDensityFSAvg(Wave.topWavelengthIndex : Wave.BSR - 1));
            AVO.density2 = mean(data.BulkDensityFSAvg(Wave.BSR : Wave.BGHSZ));
            AVO.density3 = mean(data.BulkDensityFSAvg(Wave.BGHSZ + 1 : Wave.bottomWavelengthIndex));
            
            data.BulkDensityFSAvg(Wave.topLimit : Wave.BSR - 1) = AVO.density1;
            if ~originalResolutionFlag
                data.BulkDensityFSAvg(Wave.BSR : Wave.BGHSZ) = AVO.density2;
            end
            data.BulkDensityFSAvg(Wave.BGHSZ + 1 : end) = AVO.density3;
            
            
            % Averaged VS            
            AVO.vs1 = mean(data.VSFSAvg(Wave.topWavelengthIndex : Wave.BSR - 1));
            AVO.vs2 = mean(data.VSFSAvg(Wave.BSR : Wave.BGHSZ));
            AVO.vs3 = mean(data.VSFSAvg(Wave.BGHSZ + 1 : Wave.bottomWavelengthIndex));
            
            data.VSFSAvg(Wave.topLimit : Wave.BSR - 1) = AVO.vs1;
            if ~originalResolutionFlag
                data.VSFSAvg(Wave.BSR : Wave.BGHSZ) = AVO.vs2;
            end
            data.VSFSAvg(Wave.BGHSZ + 1 : end) = AVO.vs3;
        end
        function [ impedance ] = CalcImpedance( ~ , bulkDensityFSAvg , VPFSAvg )
            impedance = bulkDensityFSAvg .* VPFSAvg;
        end
        function [ reflectionCoefficient ] = CalcReflectionCoefficient( ~ , impedance )
            I1 = impedance(1 : end - 1);
            I2 = impedance(2 : end);
            R = (I2 - I1) ...
             ./ (I2 + I1);
            
            reflectionCoefficient = [R; 0];
        end
        function [ dt ] = CalcConvolutionDt( ~ , VPFS , VPFSAvg , originalResolutionFlag )
            if originalResolutionFlag
                dt = 1 ./ VPFS;
            else
                dt = 1 ./ VPFSAvg;
            end
        end
        function [ timeSeries ] = CalcTimeSeries( ~ , convolutionDt )
            convolutionDt(isnan(convolutionDt)) = 0;
            timeSeries = cumsum(convolutionDt);
        end
        function [ F , T , C ] = CalcSeismogram( ~ , rickerFrequency , timeSeries , reflectionCoefficient )
            
            reflectionCoefficient(isnan(reflectionCoefficient)) = 0;
            
            
            f = (1 - 2 * pi ^ 2 * rickerFrequency ^ 2 .* timeSeries .^ 2) ...
                .* exp(-1 * pi ^ 2 * rickerFrequency ^ 2 .* timeSeries .^ 2);
            F = flipud(f);
            F(end + 1 : end + length(timeSeries)) = f;
            T = flipud(timeSeries .* -1);
            T(end + 1 : end + length(timeSeries)) = timeSeries;
            C = filter(reflectionCoefficient, 1, F);
            C = C(end / 2 + 1 : end);
        end
        
        
%         function [  ] = ( obj )
%             
%         end
        
        
        
        
        %%% Plotting methods
        
        function PlotPhaseSaturations( obj , Wave )
            
            figure1 = figure();
            
            quantity = [12 26 40];
            n = numel(quantity);
            
            
            
            for i = 1:n
                iQuantity = quantity(i);
                
                subplot(1, 3, i)
                hold on
                plot(obj.SaturationLF.Hydrate(:, iQuantity), Wave.depth, ...
                    'Color', [0 .5 0], ...
                    'LineWidth', 2.5);
                plot(obj.SaturationLF.Gas(:, iQuantity), Wave.depth, ...
                    'Color', [1 0 0], ...
                    'LineWidth', 2.5);
                axis([0 0.31 420 520])
                xlabel('Fluid saturations')
                ylabel('Depth (mbsf)')
                if i == 1
                    legend('Hydrate saturation', 'Gas saturation')
                end
                set(gca,'YDir','Reverse')
                
                switch i
                    case 1
                        titleString = 'a';
                    case 2
                        titleString = 'b';
                    case 3
                        titleString = 'c';
                end
                title(titleString)
            end
            
        end
        function PlotParameterSensitivity( obj , Wave , WaveBase )
            quantity = obj.selectedQuantities;
            colorStream = jet(numel(Wave.seismogram));
            
            %%%%% Figure 1
            figure1 = figure();
            
            %%% VP
            axis1 = subplot(2, 2, 4);
            hold on
            % Base case
            plot(Wave.depth, WaveBase.VPFS, 'k:', 'linewidth', 2.5)
            % Quantity cases
            for iQuantity = quantity
                plot(Wave.depth, Wave.VPFS{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
            end
            xlabel('Depth (mbsf)')
            ylabel('Compressional wave velocity (m/s)')
            axis([450 510 600 2400])
            title('d')

            %%% VS
            axis2 = subplot(2, 2, 3);
            hold on
            % Base case
            plot(Wave.depth, WaveBase.VSFS, 'k:', 'linewidth', 2.5)
            % Quantity cases
            for iQuantity = quantity
                plot(Wave.depth, Wave.VSFS{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
            end
            xlabel('Depth (mbsf)')
            ylabel('Shear wave velocity (m/s)')
            axis([450 510 400 600])
            legend( '0 g/dm^3' , '6 g/dm^3' , '15 g/dm^3' , '23 g/dm^3' , '32 g/dm^3' , '40 g/dm^3' )
            title('c')


            %%% Density
            axis3 = subplot(2, 2, 2);
            hold on
            % Base case
            plot(Wave.depth, WaveBase.bulkDensityFS ./ 1000, 'k:', 'linewidth', 2.5)
            % Quantity cases
            for iQuantity = quantity
                plot(Wave.depth, Wave.bulkDensityFS{iQuantity} ./ 1000, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
            end
            xlabel('Depth (mbsf)')
            ylabel('Bulk density (g/cm^3)')
            axis([450 510 1.5 1.8])
            title('b')

            %%% Bulk modulus
            axis4 = subplot(2, 2, 1);
            hold on
            % Base case
            plot(Wave.depth, WaveBase.bulkKFS, 'k:', 'linewidth', 2.5)
            % Quantity cases
            for iQuantity = quantity
                plot(Wave.depth, Wave.bulkKFS{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
            end
            xlabel('Depth (mbsf)')
            ylabel('Bulk modulus (Pa)')
            axis([450 510 7e8 7e9])
            title('a')
            
            axis1.Position = [.61 .09 .35 .37];
            axis2.Position = [.11 .09 .35 .37];
            axis3.Position = [.61 .58 .35 .37];
            axis4.Position = [.11 .58 .35 .37];
            
            figure1.Position(3) = 416;
            set(findall(figure1,'-property','FontSize'),'FontSize',8)
            set(findall(figure1,'-property','FontName'),'FontName','Arial')
            
            
            
            
            
            %%%%% Figure 2
            figure2 = figure();
            
            %%% Depth series seismogram
            axis5 = subplot(2, 1, 1);
            hold on
            % 3 phase bulk equilibrium depth line
            BEQL3P = 481; % mbsf
            plot([BEQL3P, BEQL3P], [-1, 1], '--' , 'Color' , [.4 .4 .4] , 'linewidth' , 1.5 );
            
            figureCellArray = cell(numel(quantity), 1);
            figureNumber = 0;
            
            for iQuantity = quantity
                figureNumber = figureNumber + 1;
                figureCellArray{figureNumber} = plot(Wave.depth, Wave.seismogram{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
            end
            axis([380 580 -.35 .2])
            xlabel('Depth (mbsf)')
            ylabel('Amplitude')
            title('a) Depth Series')
            legend( [figureCellArray{:}], '6 g/dm^3' , '15 g/dm^3' , '23 g/dm^3' , '32 g/dm^3' , '40 g/dm^3' )
            
            
            %%% Time series seismogram
            axis6 = subplot(2, 1, 2);
            hold on
            for iQuantity = quantity
                plot(Wave.time{iQuantity} - Wave.time{iQuantity}(Wave.depth == 380) , Wave.seismogram{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5)
            end       
            axis([0 Wave.time{iQuantity}(Wave.depth == 580) - Wave.time{iQuantity}(Wave.depth == 380) -.35 .2])
            xlabel('Time (seconds)')
            ylabel('Amplitude')
            title('b) Time Series')
            legend( '6 g/dm^3' , '15 g/dm^3' , '23 g/dm^3' , '32 g/dm^3' , '40 g/dm^3' )
            
            axis5.Position = [.13 .56 .79 .39];
            axis6.Position = [.13 .07 .79 .39];
            
            figure2.Position(3) = 416;
            set(findall(figure2,'-property','FontSize'),'FontSize',8)
            set(findall(figure2,'-property','FontName'),'FontName','Arial')
        end
        function PlotBackgroundProperties( obj , dataBase )
            
            BSR_top_estimate = 456; % mbsf
            BSR_bottom_estimate = 472; % mbsf
            three_phase_bulk_equilibrium = 481; % mbsf

            BSR_top_line = [BSR_top_estimate BSR_top_estimate];
            BSR_bottom_line = [BSR_bottom_estimate BSR_bottom_estimate];
            bulk_EQL_line = [three_phase_bulk_equilibrium three_phase_bulk_equilibrium];

            plotLineWidth = 1.5;
            markerLineWidth = 1.25;






            figure1 = figure();

            figure1.Position(3) = 800;
            set(findall(figure1,'-property','FontSize'),'FontSize',8)
            set(findall(figure1,'-property','FontName'),'FontName','Arial')


            %%% Gamma ray
            leftLimit = 12;
            rightLimit = 110;    
            xLine = [leftLimit rightLimit];

            axis1 = subplot(1,7,1);

            rectangle('Position', [leftLimit, BSR_top_estimate, rightLimit - leftLimit, BSR_bottom_estimate - BSR_top_estimate], ...
                        'FaceColor', [0.8, 0.8, 0.8], ...
                        'EdgeColor', 'k');
            hold on

            xText = leftLimit + 0.02*(rightLimit - leftLimit);
            text(xText, (BSR_top_estimate + BSR_bottom_estimate)/2, 'BSR')
            text(xText, bulk_EQL_line(1) + (BSR_bottom_estimate - BSR_top_estimate)/2, '3P EQL')



            plot(dataBase.GammaRay, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
            xlabel('Gamma ray (gAPI)','FontSize',8)

            axis([leftLimit rightLimit 160 500])
            title('a')    
            hold on
        %     plot( x_line , BSR_top_line , 'k--' , 'linewidth' , line_width )
        %     hold on
        %     plot( x_line , BSR_bottom_line , 'k--' , 'linewidth' , line_width )
        %     hold on
            plot( xLine , bulk_EQL_line , 'k--' , 'linewidth' , markerLineWidth )    

            axis1.FontSize = 8;
            axis1.XTick = [30 65 100];




            % Caliper
            leftLimit = 4;
            rightLimit = 10;    
            xLine = [leftLimit rightLimit];

            bit_diameter = 10.125; % inches
            bit_radius = bit_diameter ./ 2;
            CALI.Radius = obj.CALI.Diameter ./ 2;

            axis2 = subplot(1,7,2);

            rectangle('Position', [leftLimit, BSR_top_estimate, rightLimit - leftLimit, BSR_bottom_estimate - BSR_top_estimate], ...
                        'FaceColor', [0.8, 0.8, 0.8], ...
                        'EdgeColor', 'k');
            hold on
%             plot(dataBase.GammaRay, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
            plot( CALI.Radius , obj.CALI.Depth , 'k' , 'linewidth' , .5 )
            hold on
            plot( [bit_radius bit_radius] , [160 500] , 'k-.' , 'linewidth' , 1 )



        %     plot( -CALI.Radius , CALI.Depth , 'k' , 'linewidth' , plot_width )

            set(gca,'YDir','Reverse')
            xlabel('Hole radius (in)')
            axis([leftLimit rightLimit 160 500])
            title('b')    
            hold on
        %     plot( x_line , BSR_top_line , 'k--' , 'linewidth' , line_width )
        %     hold on
        %     plot( x_line , BSR_bottom_line , 'k--' , 'linewidth' , line_width )
        %     hold on
            plot( xLine , bulk_EQL_line , 'k--' , 'linewidth' , markerLineWidth )    





            axis2.FontSize = 8;
            axis2.YTickLabel = [];
            axis2.XTick = [5 7 9];








            % Resisitivity
            leftLimit = 0.5;
            rightLimit = 1.5;
            xLine = [leftLimit rightLimit];

            axis3 = subplot(1,7,3);

            rectangle('Position', [leftLimit, BSR_top_estimate, rightLimit - leftLimit, BSR_bottom_estimate - BSR_top_estimate], ...
                        'FaceColor', [0.8, 0.8, 0.8], ...
                        'EdgeColor', 'k');
            hold on

            plot(dataBase.Resistivity, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
%             plot( Data.log(:,21) , Data.log(:,1) , 'k' , 'linewidth' , plotLineWidth )
            set(gca,'YDir','Reverse')
            xlabel('Resistivity (ohm-m)')
            axis([leftLimit rightLimit 160 500])
            title('c')
            hold on
        %     plot( x_line , BSR_top_line , 'k--' , 'linewidth' , line_width )
        %     hold on
        %     plot( x_line , BSR_bottom_line , 'k--' , 'linewidth' , line_width )
        %     hold on
            plot( xLine , bulk_EQL_line , 'k--' , 'linewidth' , markerLineWidth )

            axis3.FontSize = 8;
            axis3.YTickLabel = [];
            axis3.XTick = [0.6 1 1.4];






            % Calculated porosity
            leftLimit = 0.4;
            rightLimit = 0.8;
            xLine = [leftLimit rightLimit];

            axis4 = subplot(1,7,4);

            rectangle('Position', [leftLimit, BSR_top_estimate, rightLimit - leftLimit, BSR_bottom_estimate - BSR_top_estimate], ...
                        'FaceColor', [0.8, 0.8, 0.8], ...
                        'EdgeColor', 'k');
            hold on    
            plot(dataBase.Porosity, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
%             plot( Data.log(:,3) , Data.log(:,1) , 'k' , 'linewidth' , plotLineWidth )
            set(gca,'YDir','Reverse')
            xlabel('Porosity')
            axis([leftLimit rightLimit 160 500])
            title('d')  
            hold on
        %     plot( x_line , BSR_top_line , 'k--' , 'linewidth' , line_width )
        %     hold on
        %     plot( x_line , BSR_bottom_line , 'k--' , 'linewidth' , line_width )
        %     hold on
            plot( xLine , bulk_EQL_line , 'k--' , 'linewidth' , markerLineWidth )


            axis4.FontSize = 8;
            axis4.YTickLabel = [];
            axis4.XTick = [0.45 0.6 0.75];
            
            
            
            
            % Bulk density
            leftLimit = 1.3;
            rightLimit = 1.9;    
            xLine = [leftLimit rightLimit];

            axis5 = subplot(1,7,5);


            rectangle('Position', [leftLimit, BSR_top_estimate, rightLimit - leftLimit, BSR_bottom_estimate - BSR_top_estimate], ...
                        'FaceColor', [0.8, 0.8, 0.8], ...
                        'EdgeColor', 'k');
            hold on

            xText = leftLimit + 0.02*(rightLimit - leftLimit);
            text(xText, (BSR_top_estimate + BSR_bottom_estimate)/2, 'BSR')
            text(xText, bulk_EQL_line(1) + 0.6*(BSR_bottom_estimate - BSR_top_estimate), '3P EQL')    


            plot(dataBase.BulkDensity ./ 1000, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
%             plot( Data.log(:,6) , Data.log(:,1) , 'k' , 'linewidth' , plotLineWidth )
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
            xlabel('Bulk density (g/cm^3)')
            axis([leftLimit rightLimit 160 500])
            title('e')
            hold on
        %     plot( x_line , BSR_top_line , 'k--' , 'linewidth' , line_width )
        %     hold on
        %     plot( x_line , BSR_bottom_line , 'k--' , 'linewidth' , line_width )
        %     hold on
            plot( xLine , bulk_EQL_line , 'k--' , 'linewidth' , markerLineWidth )    

            axis5.FontSize = 8;


            % VP
            leftLimit = 1400;
            rightLimit = 2000;    
            xLine = [leftLimit rightLimit];

            axis6 = subplot(1,7,6);

            rectangle('Position', [leftLimit, BSR_top_estimate, rightLimit - leftLimit, BSR_bottom_estimate - BSR_top_estimate], ...
                        'FaceColor', [0.8, 0.8, 0.8], ...
                        'EdgeColor', 'k');
            hold on    
            plot(dataBase.VP, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
%             plot( Data.log(:,4) , Data.log(:,1) , 'k' , 'linewidth' , plotLineWidth )
            set(gca,'YDir','Reverse')
            xlabel('P-velocity (m/s)')
            axis([leftLimit rightLimit 160 500])
            title('f')
            hold on
        %     plot( x_line , BSR_top_line , 'k--' , 'linewidth' , line_width )
        %     hold on
        %     plot( x_line , BSR_bottom_line , 'k--' , 'linewidth' , line_width )
        %     hold on
            plot( xLine , bulk_EQL_line , 'k--' , 'linewidth' , markerLineWidth )    

            axis6.FontSize = 8;
            axis6.YTickLabel = [];
            axis6.XTick = [1500 1700 1900];



            % VS
            leftLimit = 100;
            rightLimit = 900;    
            xLine = [leftLimit rightLimit];

            axis7 = subplot(1,7,7);

            rectangle('Position', [leftLimit, BSR_top_estimate, rightLimit - leftLimit, BSR_bottom_estimate - BSR_top_estimate], ...
                        'FaceColor', [0.8, 0.8, 0.8], ...
                        'EdgeColor', 'k');
            hold on
            plot(dataBase.VS, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
%             plot( Data.log(:,5) , Data.log(:,1) , 'k' , 'linewidth' , plotLineWidth )
            set(gca,'YDir','Reverse')
            xlabel('S-velocity (m/s)')
            axis([leftLimit rightLimit 160 500])
            title('g')
            hold on
        %     plot( x_line , BSR_top_line , 'k--' , 'linewidth' , line_width )
        %     hold on
        %     plot( x_line , BSR_bottom_line , 'k--' , 'linewidth' , line_width )
        %     hold on
            plot( xLine , bulk_EQL_line , 'k--' , 'linewidth' , markerLineWidth )

            axis7.FontSize = 8;
            axis7.YTickLabel = [];
            axis7.XTick = [250 500 750];


            interval = 0.13;
            gap = 0.01;
            xx = 0.06;
            axis1.Position(1) = xx;
            axis1.Position(3) = interval - gap;
            xx = xx + interval;
            axis2.Position(1) = xx;
            axis2.Position(3) = interval - gap;
            xx = xx + interval;
            axis3.Position(1) = xx;
            axis3.Position(3) = interval - gap;
            xx = xx + interval;
            axis4.Position(1) = xx;
            axis4.Position(3) = interval - gap;            
            xx = xx + interval;
            axis5.Position(1) = xx;
            axis5.Position(3) = interval - gap;
            xx = xx + interval;
            axis6.Position(1) = xx;
            axis6.Position(3) = interval - gap;
            xx = xx + interval;
            axis7.Position(1) = xx;
            axis7.Position(3) = interval - gap;
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























