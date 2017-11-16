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
        function [ Wave , data ] = RunSeismicAnalysisRoutine( obj , caseString )
            %%% Make main table
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
            
            % Density in kg/m^3
            data.BulkDensity = nan(nDepth, 1);
            data.BulkDensity(obj.LDT.Fixed_Depth_TOP : obj.LDT.Fixed_Depth_BOTTOM) = obj.LDT.Fixed_RHOB .* 1000;

            
            data.HydrateK = obj.CalcHydrateK(data.Pressure, data.Temperature);
            
            
            testFlag = true;
%             testFlag = false;
            if testFlag
                data.Porosity = obj.CalcPorosity(data.Resistivity, 1, data.Temperature);
            end
            
            
            
            
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
                    
%                     data.HydrateK(:) = mean(data.HydrateK(upperBoundIndex : lowerBoundIndex));
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
            
            
            for iQuantity = 1:nQuantity
                quantity = obj.quantityArray(iQuantity);
                
                [ seismogram , timeSeries , thickness , parameterSensitivity ] ...
                    = obj.RunRockPhysicsRoutine(data, Wave, iQuantity, caseString, testFlag);
                
                
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
                RunRockPhysicsRoutine( obj , data , Wave , iQuantity , caseString , testFlag )
            
            seismogram = [];
            timeSeries = [];
            thickness = [];
            parameterSensitivity = [];
            
            
            
            
            
            
            data.Sh = obj.SaturationLF.Hydrate(:, iQuantity);
            data.Sg = obj.SaturationLF.Gas(:, iQuantity);
            data.Sw = 1 - data.Sh - data.Sg;
            % Gets indices of depths at BSR and BGHSZ
            Wave.BSR = find(data.Depth == obj.graphTop + obj.SaturationLF.Top3P(iQuantity));
            Wave.BGHSZ = find(data.Depth == obj.graphTop + obj.SaturationLF.Bottom3P(iQuantity));
            Wave.rickerFrequency = 30;
            
            if ~testFlag
                data.Porosity = obj.CalcPorosity(data.Resistivity, data.Sw);
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
            R = (I1 - I2) ...
             ./ (I1 + I2);
            
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
        %
        function PlotPhaseSaturations( obj , Wave )
            
            figObj = figure;
            
%             quantity = [6 23 40];
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
%                 axis([0 0.2 420 520])
                axis([0 0.31 420 520])
                xlabel('Fluid saturations')
                ylabel('Depth (mbsf)')
                if i == 1
                    legend('Hydrate saturation', 'Gas saturation')
                end
                set(gca,'YDir','Reverse')
                
%                 titleString = strcat('Methane Quantity = ', num2str(iQuantity), ' (g/dm^3 of pore volume)');
                switch i
                    case 1
                        titleString = 'a)';
                    case 2
                        titleString = 'b)';
                    case 3
                        titleString = 'c)';
                end
                title(titleString)
            end
            
        end
        function PlotParameterSensitivity( obj , Wave )
            quantity = obj.selectedQuantities;
            
            figure
            
            
            colorStream = colormap(jet(numel(Wave.seismogram)));
            
            %%% VP
            axis1 = subplot(2, 2, 4);
            hold on
            % Hard coded Sw = 1 case for the 2nd forloop
%             plot( Wave.depth , %%%Data.VP(:,2) , 'k:' , 'linewidth' , 2.5 )
            
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
            % Hard coded Sw = 1 case for the 2nd forloop
%             plot( Data.log(:,1) , Data.VS(:,2) , 'k:' , 'linewidth' , 2.5 )
            for iQuantity = quantity
                plot(Wave.depth, Wave.VSFS{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
%                 plot( Data.log(:,1) , Data.VS(:,iQuantity) , 'Color' , colorStream(iQuantity,:)' , 'linewidth' , 2.5 );
            end
            xlabel('Depth (mbsf)')
            ylabel('Shear wave velocity (m/s)')
            axis([450 510 400 600])
            legend( '0 g/dm^3' , '6 g/dm^3' , '15 g/dm^3' , '23 g/dm^3' , '32 g/dm^3' , '40 g/dm^3' )
            title('c')


            %%% Density
            axis3 = subplot(2, 2, 2);
            hold on
%             % Hard coded Sw = 1 case for the 2nd forloop
%             plot( Data.log(:,1) , Data.Density(:,2) , 'k:' , 'linewidth' , 2.5 )
            for iQuantity = quantity
                plot(Wave.depth, Wave.bulkDensityFS{iQuantity} ./ 1000, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
%                 plot( Data.log(:,1) , Data.Density(:,iQuantity) , 'Color' , colorStream(iQuantity,:)' , 'linewidth' , 2.5 );
            end
            xlabel('Depth (mbsf)')
            ylabel('Bulk density (g/cm^3)')
            axis([450 510 1.5 1.8])
            title('b')

            % K Bulk Modulus
            % figure
            % plot( Data.log(:,1) , Data.log(:,8) , 'k.' , 'linewidth' , 2.5 )
            % hold on
            axis4 = subplot(2, 2, 1);
            hold on
            % Hard coded Sw = 1 case for the 2nd forloop
%             plot( Data.log(:,1) , Data.K_bulk(:,2) , 'k:' , 'linewidth' , 2.5 )
%             colorStream = colormap(jet(size(Wave.seismogram_data,2)));
            for iQuantity = quantity
                plot(Wave.depth, Wave.bulkKFS{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
%                 plot( Data.log(:,1) , Data.K_bulk(:,iQuantity) , 'Color' , colorStream(iQuantity,:)' , 'linewidth' , 2.5 );
            end
            xlabel('Depth (mbsf)')
            ylabel('Bulk modulus (Pa)')
            axis([450 510 7e8 7e9])
            title('a')



            axis1.Position = [.61 .09 .35 .37];
            axis2.Position = [.11 .09 .35 .37];
            axis3.Position = [.61 .58 .35 .37];
            axis4.Position = [.11 .58 .35 .37];

            % c3 = colorbar('Limits',[6 Max_Concentration],...
            %     'Ticks',linspace(6,40,5),...
            %     'TickLabels',linspace(6,40,5));
            % caxis([6 Max_Concentration])
            % c3.Label.String = 'Methane quantity (g/dm^3 of pore volume)';
            % c3.Label.FontSize = 11;



            figure1 = gcf;
            figure1.Position(3) = 416;
            set(findall(figure1,'-property','FontSize'),'FontSize',8)
            set(findall(figure1,'-property','FontName'),'FontName','Arial')
        end
        %}
        
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























