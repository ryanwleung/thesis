classdef DCSeismicAnalysisBR < DCBlakeRidge
    properties
        graphTop
        graphBottom
        saturationBottom
        
        depthArrayFull
        
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
        
        seaVP = 1485; % m/s
        
        selectedQuantities = [6 15 23 32 40];
        
        aHelgerudK = -1.09e-2; % GPa/C deg
        bHelgerudK = 3.8e-3; % GPa/MPa
        cHelgerudK = 8.39; % GPa
        
        phiReferenceRo = 0.24 % ohm-m (.223)
        phiReferenceTemp = 18; % C deg
        phiA = 0.9;
        phiM = 2.7;
        phiN = 1.9386;
        
%         gasEOS = 'IdealGasLaw';
        gasEOS = 'Batzle';
        
        rickerFrequency = 30;
        
        searchPeakTopIndex = 3200;
        searchPeakBottomIndex = 3300;
        
        
        axisMaxAmplitude = 0.075;
        axisMinAmplitude = -0.1;
    end
    methods
        %%% Constructor
        function [ obj ] = DCSeismicAnalysisBR()
            obj@DCBlakeRidge()
            
            obj.graphTop = obj.liuTop - obj.seafloorDepth; % mbsf
            obj.graphBottom = obj.liuBottom - obj.seafloorDepth; % mbsf
            obj.saturationBottom = obj.maxDepth; % mbsf

            obj.depthArray = (obj.minDepth : 1 : obj.maxDepth)';
            obj.depthArrayFull = [(1 : 1 : obj.seafloorDepth)'; obj.depthArray + obj.seafloorDepth];
            
            obj.quantityArray = (1 : 1 : 40)';
            
            obj.SaturationLF = obj.LoadPhaseBehaviorBlakeRidge();
        end
        
        %%% Main methods
        function [ Wave , data , WaveBase , dataBase , WaveDickens ] = RunSeismicAnalysisRoutine( obj , caseString )
            %%% Instantiate optional base (Sw = 1) calculation variables
            WaveBase = struct();
            dataBase = table();
            
            %%% Instantiate Wave structs
            Wave = struct();
            WaveDickens = struct();
            
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
            
            data.Porosity = obj.CalcPorosity(data.Resistivity, 1, data.Temperature);
            
            %%% Create table for background property figure with porosity
            %%% at Sw = 1
            dataBase = data;
            dataBase.Porosity = obj.CalcPorosity(data.Resistivity, 1, data.Temperature);
            %%% end
            
            
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
            
            Wave.bulkKFS = cell(nQuantity, 1);
            Wave.bulkDensityFS = cell(nQuantity, 1);
            Wave.VSFS = cell(nQuantity, 1);
            Wave.VPFS = cell(nQuantity, 1);
            
            
            %%% Run base case for parameter sensitivity
            [ ~ , ~ , ~ , WaveBase ] ...
                    = obj.RunRockPhysicsRoutine(data, Wave, 0, caseString);
            
            %%% For loop through all methane quantities
            for iQuantity = 1:nQuantity                
                %%% Get seismogram and related outputs from routine
                [ seismogram , timeSeries , thickness , parameterSensitivity ] ...
                    = obj.RunRockPhysicsRoutine(data, Wave, iQuantity, caseString);
                
                
                Wave.seismogram{iQuantity} = seismogram;
                Wave.time{iQuantity} = timeSeries;
                Wave.thickness3P{iQuantity} = thickness;
                
                %%% Find the leading and trailing peak amplitude
                [Wave.peak{iQuantity, 1}, Wave.peak{iQuantity, 2}] = obj.GetPeakAmplitudes(seismogram);
                
                switch caseString
                    case 'ParameterSensitivity'
                        Wave.bulkKFS{iQuantity} = parameterSensitivity.bulkKFS;
                        Wave.bulkDensityFS{iQuantity} = parameterSensitivity.bulkDensityFS;
                        Wave.VSFS{iQuantity} = parameterSensitivity.VSFS;
                        Wave.VPFS{iQuantity} = parameterSensitivity.VPFS;
                    case 'OriginalResolution'
                        Wave.VPFS{iQuantity} = parameterSensitivity.VPFS;
                end
            end
            
            %%% Run Dickens case
            if ~strcmp('ParameterSensitivity', caseString)
                [ seismogram , timeSeries , ~ , ~ ] ...
                        = obj.RunRockPhysicsRoutine(data, Wave, -1, caseString);
                WaveDickens.seismogram = seismogram;
                WaveDickens.time = timeSeries;
            end
        end
        function [ seismogram , timeSeries , thickness , parameterSensitivity ] = ...
                RunRockPhysicsRoutine( obj , data , Wave , iQuantity , caseString )
            
            seismogram = [];
            timeSeries = [];
            thickness = -100;
            parameterSensitivity = [];
            
            
            Wave.BSR = [];
            Wave.BGHSZ = [];
            switch iQuantity
                case 0
                    data.Sh = zeros(height(data), 1);
                    data.Sg = zeros(height(data), 1);
                    data.Sw = ones(height(data), 1);
                case -1
                    data.Sh = obj.Dickens.Hydrate;
                    data.Sg = obj.Dickens.Gas;
                    data.Sw = 1 - data.Sh - data.Sg;
                    
                    % Sets indices of depths at BSR and BGHSZ
%                     Wave.BSR = 354;
%                     Wave.BGHSZ = 376;
                    Wave.BSR = 469;
                    Wave.BGHSZ = 490;
                otherwise
                    data.Sh = obj.SaturationLF.Hydrate(:, iQuantity);
                    data.Sg = obj.SaturationLF.Gas(:, iQuantity);
                    data.Sw = 1 - data.Sh - data.Sg;

                    % Gets indices of depths at BSR and BGHSZ
                    Wave.BSR = find(data.Depth == obj.graphTop + obj.SaturationLF.Top3P(iQuantity));
                    Wave.BGHSZ = find(data.Depth == obj.graphTop + obj.SaturationLF.Bottom3P(iQuantity));
            end
            
            % Gas bulk modulus in Pa
            data.GasK = obj.CalcGasK(data.Pressure, data.Temperature);
            
            switch caseString
                case 'ParameterSensitivity'
                    upperBoundIndex = 415;
                    lowerBoundIndex = 630;
                    data.Porosity(:) = mean(data.Porosity(upperBoundIndex : lowerBoundIndex));
                    data.GasK(:) = mean(data.GasK(upperBoundIndex : lowerBoundIndex));
            end
            
            
            % Bulk modulus from VP and VS in Pa
            data.BulkK = obj.CalcBulkK(data.BulkDensity, data.VP, data.VS);
            data.FluidK = obj.CalcFluidK(data.Sw, obj.waterK, data.Sg, data.GasK, data.Sh, data.HydrateK);
            data.GrainK = obj.CalcGrainK(data.GammaRay);
            data.FrameK = obj.CalcFrameK(data.Porosity, data.GrainK);
            data.ShearG = obj.CalcShearG(data.BulkDensity, data.VS);

            %%% Fluid substitution
            data.BulkKFS = obj.CalcBulkKFluidSubstituted(data.Porosity, data.FluidK, data.GrainK, data.FrameK);
            data.BulkDensityFS = obj.CalcBulkDensityFluidSubstituted(data.GammaRay, data.Porosity, data.Sw, data.Sg, data.Sh);
            data.VSFS = obj.CalcVSFS(data.ShearG, data.BulkDensityFS);
            data.VPFS = obj.CalcVPFS(data.BulkKFS, data.BulkDensityFS, data.VSFS);
            
            %%% Fill in VP and bulk density values in the fluid substituted
            %%% columns for full reflectance series, using best fit through
            %%% 0-200 mbsf log points
            data.VPFS = obj.FillInSeafloorToSaturationTop(data.VPFS);
            data.BulkDensityFS = obj.FillInSeafloorToSaturationTop(data.BulkDensityFS);
            
            
            
            Wave = obj.FindIntervalOfOneSeismicWavelength(data.VPFS, Wave);
            
            switch caseString
                case 'OriginalResolution'
                    originalResolutionFlag = true;
                otherwise
                    % it appears that I
                    % always use the original resolution (no averaging of
                    % values within the 3P zone
                    originalResolutionFlag = true;
            end
            [ data , Wave ] = obj.CalcAverageProperties(data, Wave, originalResolutionFlag);
            
            %%% New table to hold values from sea level
            dataFull = table();
            seaValues = obj.GenerateSeaValues();
            
            dataFull.Depth = obj.depthArrayFull;
            dataFull.VPFS = [seaValues.VP; data.VPFS];
            dataFull.VPFSAvg = [seaValues.VP; data.VPFSAvg];
            dataFull.BulkDensityFSAvg = [seaValues.bulkDensity; data.BulkDensityFSAvg];
            
            
            dataFull.Impedance = obj.CalcImpedance(dataFull.BulkDensityFSAvg, dataFull.VPFSAvg);
            dataFull.ReflectionCoefficient = obj.CalcReflectionCoefficient(dataFull.Impedance);
            dataFull.ConvolutionDt = obj.CalcConvolutionDt(dataFull.VPFS, dataFull.VPFSAvg, originalResolutionFlag);
            dataFull.TimeSeries = obj.CalcTimeSeries(dataFull.ConvolutionDt);
            
            [ ~ , ~ , seismogram ] = obj.CalcSeismogram(dataFull.TimeSeries, dataFull.ReflectionCoefficient);
            timeSeries = dataFull.TimeSeries;
            thickness = data.Depth(Wave.BGHSZ) - data.Depth(Wave.BSR);
            if isempty(thickness)
                thickness = -100;
            end
            switch caseString
                case 'ParameterSensitivity'
                    parameterSensitivity.bulkKFS = data.BulkKFS;
                    parameterSensitivity.bulkDensityFS = data.BulkDensityFS;
                    parameterSensitivity.VSFS = data.VSFS;
                    parameterSensitivity.VPFS = data.VPFS;
                case 'OriginalResolution'
                    parameterSensitivity.VPFS = data.VPFS;
            end
        end
        
        %%% Rock physics calculation methods
        function [ temperature ] = CalcTemperature995( obj )
            depth = obj.depthArray ./ 1000; % convert to km
            temperature = obj.seafloorTemperature995 + depth .* obj.temperatureGradient995; % in C deg
        end
        function [ PorosityCorrected ] = CalcPorosity( obj , Rt , Sw , temperature )
            referenceRo = obj.phiReferenceRo;
            referenceTemp = obj.phiReferenceTemp;
            
            a = obj.phiA;
            m = obj.phiM;
            n = obj.phiN;
            
            Rw = referenceRo .* (referenceTemp + 21.5) ./ (temperature + 21.5);
            Porosity = (a .* Rw ./ Rt ...
                            ./ (Sw .^ n)) ...
                            .^ (1 ./ m);
%             sum(Sw ~= 1)
            %%% Porosity correction
            Ro = 0.8495 + 2.986e-4 .* obj.depthArray;
            
            sh = 1 - (Ro ./ Rt) .^ (1 / n);
%             sum((sh < 0))
            sh(sh < 0) = 0;
            
            PorosityCorrected = Porosity ./ ...
                                (1 - sh);
        end
        function [ gasK ] = CalcGasK( obj , pressure , temperature )
            
            
            switch obj.gasEOS
                case 'IdealGasLaw'
                    gasK = pressure;
                case 'Batzle'
                    pr1978Obj = BCPR1978Util();
                    
                    reducedPressure = pressure ./ pr1978Obj.pcCH4;
                    reducedTemperature = (temperature + 273.15) ./ pr1978Obj.tcCH4;
                    
                    gasK = obj.CalcKsBatzle(reducedPressure, reducedTemperature, pressure);

                    
                    
            end
        end
        function [ Ks ] = CalcKsBatzle( ~ , pr , tr , pressure )
            
            E = 0.109 .* (3.85 - tr) .^ 2 ...
                .* exp( -(0.45 + 8 .* (0.56 - 1 ./ tr) .^ 2) ...
                        .* pr .^ 1.2 ...
                        ./ tr ...
                      );
            
            Z = pr .* (0.03 + 0.00527 .* (3.5 - tr) .^ 3) ...
                + (0.642 .* tr ...
                    - 0.007 .* tr .^ 4 ...
                    - 0.52) ...
                + E;
            
            gamma0 = 0.85 ...
                        + 5.6 ./ (pr + 2) ...
                        + 27.1 ./ (pr + 3.5) .^ 2 ...
                        - 8.7 .* exp(-0.65 .* (pr + 1));
            
            dZdPr = (0.03 + 0.00527 .* (3.5 - tr) .^ 3) ...
                    + (0.109 .* (3.85 - tr) .^ 2) ...
                        .* (1.2 .* pr .^ 0.2) ...
                            .* (-(0.45 + 8 .* (0.56 - 1 ./ tr) .^ 2) ./ tr) ...
                        .* exp( -(0.45 + 8 .* (0.56 - 1 ./ tr) .^ 2) ...
                                .* pr .^ 1.2 ...
                                ./ tr ...
                              );
            
            Ks = pressure ...
                    .* gamma0 ...
                    ./ (1 - pr ./ Z .* dZdPr);
            
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
            % Fluid modulus using isostress average (volume-weighted harmonic mean) 
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
        function [ Wave ] = FindIntervalOfOneSeismicWavelength( obj , VPFS , Wave )
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
            
            frequency = obj.rickerFrequency;
            
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
                topWavelength = topVelocity / frequency;
                
                bottomVelocity = mean(VPFS(Wave.BGHSZ + 1 : Wave.BGHSZ + iterationThickness));
                bottomWavelength = bottomVelocity / frequency;
                
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
            vp1 = mean(data.VPFSAvg(Wave.topWavelengthIndex : Wave.BSR - 1));
            vp2 = mean(data.VPFSAvg(Wave.BSR : Wave.BGHSZ));
            vp3 = mean(data.VPFSAvg(Wave.BGHSZ + 1 : Wave.bottomWavelengthIndex));
            
            data.VPFSAvg(Wave.topLimit : Wave.BSR - 1) = vp1;
            if ~originalResolutionFlag
                data.VPFSAvg(Wave.BSR : Wave.BGHSZ) = vp2;
            end
            data.VPFSAvg(Wave.BGHSZ + 1 : end) = vp3;
            
            
            % Averaged Density
            density1 = mean(data.BulkDensityFSAvg(Wave.topWavelengthIndex : Wave.BSR - 1));
            density2 = mean(data.BulkDensityFSAvg(Wave.BSR : Wave.BGHSZ));
            density3 = mean(data.BulkDensityFSAvg(Wave.BGHSZ + 1 : Wave.bottomWavelengthIndex));
            
            data.BulkDensityFSAvg(Wave.topLimit : Wave.BSR - 1) = density1;
            if ~originalResolutionFlag
                data.BulkDensityFSAvg(Wave.BSR : Wave.BGHSZ) = density2;
            end
            data.BulkDensityFSAvg(Wave.BGHSZ + 1 : end) = density3;
            
            
            % Averaged VS            
            vs1 = mean(data.VSFSAvg(Wave.topWavelengthIndex : Wave.BSR - 1));
            vs2 = mean(data.VSFSAvg(Wave.BSR : Wave.BGHSZ));
            vs3 = mean(data.VSFSAvg(Wave.BGHSZ + 1 : Wave.bottomWavelengthIndex));
            
            data.VSFSAvg(Wave.topLimit : Wave.BSR - 1) = vs1;
            if ~originalResolutionFlag
                data.VSFSAvg(Wave.BSR : Wave.BGHSZ) = vs2;
            end
            data.VSFSAvg(Wave.BGHSZ + 1 : end) = vs3;
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
            
            if any(isnan(dt))
                temp = dt;
                temp(isnan(dt)) = [];
                dt(isnan(dt)) = temp(end);
            end
        end
        function [ timeSeries ] = CalcTimeSeries( ~ , convolutionDt )
            timeSeries = cumsum(convolutionDt);
            timeSeries = timeSeries .* 2;
        end
        function [ F , T , C ] = CalcSeismogram( obj , timeSeries , reflectionCoefficient )
            logical = isnan(reflectionCoefficient);
%             if sum(logical) ~= 0
%                 sum(logical)
%             end
            reflectionCoefficient(logical) = 0;
            
            % Reduce frequency (double wavelength of wavelet) to account
            % for TWT time convolution
            frequency = obj.rickerFrequency / 2;
            
            %%% Original code from Michael's rock physics project
            f = (1 - 2 * pi ^ 2 * frequency ^ 2 .* timeSeries .^ 2) ...
                .* exp(-1 * pi ^ 2 * frequency ^ 2 .* timeSeries .^ 2);
            F = flipud(f);
            F(end + 1 : end + length(timeSeries)) = f;
            T = flipud(timeSeries .* -1);
            T(end + 1 : end + length(timeSeries)) = timeSeries;
            C = filter(reflectionCoefficient, 1, F);
            C = C(end / 2 + 1 : end);
        end
        function [ adjustedDataSeries ] = FillInSeafloorToSaturationTop( obj , dataSeries )
            topDepth = 1;
            bottomDepth = 199;
            
            x = obj.depthArray(topDepth : bottomDepth);
            y = dataSeries(topDepth : bottomDepth);
            nanLogical = isnan(y);
            polyfitDegree1 = polyfit(x(~nanLogical), y(~nanLogical), 1);
            
            adjustedDataSeries = dataSeries;
            adjustedDataSeries(nanLogical) = polyval(polyfitDegree1, x(nanLogical));
        end
        function [ seaValues ] = GenerateSeaValues( obj )
            seaValues = struct();
            
            n = obj.seafloorDepth;
            
            seaValues.VP = obj.seaVP .* ones(n, 1);
            seaValues.bulkDensity = obj.waterDensity .* ones(n, 1);
        end
        function [ leadingAmp , trailingAmp ] = GetPeakAmplitudes( obj , seismogram )
            leadingAmp = [];
            trailingAmp = [];
            
            [peaks, indices, widths, prominences] = findpeaks(seismogram);
%             findpeaks(seismogram)
            
            peaksNearBSR = peaks(indices >= obj.searchPeakTopIndex & indices <= obj.searchPeakBottomIndex);
            
            if numel(peaksNearBSR) == 2
                leadingAmp = peaksNearBSR(1);
                trailingAmp = peaksNearBSR(2);
            end
        end
        
        
        
        
        
        
        %%% Plotting methods
        function PlotPhaseSaturations( obj )
            
            figure1 = figure();
            
            quantity = [12 26 40];
            n = numel(quantity);
            
            depth = obj.depthArray + obj.seafloorDepth;
            
            for i = 1:n
                iQuantity = quantity(i);
                
                axisCell = subplot(1, 3, i);
                hold on
                plot(obj.SaturationLF.Hydrate(:, iQuantity), depth, ...
                    'Color', [0 .5 0], ...
                    'LineWidth', 2.5);
                plot(obj.SaturationLF.Gas(:, iQuantity), depth, ...
                    'Color', [1 0 0], ...
                    'LineWidth', 2.5);
                axis([0 0.31 420 + obj.seafloorDepth 520 + obj.seafloorDepth])
                xlabel('Fluid saturations')
                if i == 1
                    ylabel('Depth (mbsl)')
                else
                    axisCell.YTickLabel = [];
                end
                if i == 3
                    legend('Hydrate', 'Gas')
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
            
            depth = obj.depthArray + obj.seafloorDepth;
            
            axisMinDepth = 450 + obj.seafloorDepth;
            axisMaxDepth = 510 + obj.seafloorDepth;
            
            leftLimit = 3248;
            rightLimit = 3276;
            
            
            %%%%% Figure 1
            figure1 = figure();
            
            %%% VP
            axis1 = subplot(2, 2, 4);
            hold on
            
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, 1400, 2400)
            
            % Base case
            plot(depth, WaveBase.VPFS, 'k:', 'linewidth', 2.5)
            % Quantity cases
            for iQuantity = quantity
                plot(depth, Wave.VPFS{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
            end
            xlabel('Depth (mbsl)')
            ylabel('Compressional wave velocity (m/s)')
            axis([axisMinDepth axisMaxDepth 1400 2400])
            title('d')

            %%% VS
            axis2 = subplot(2, 2, 3);
            hold on
            
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, 530, 600)
            
            % Base case
            plot(depth, WaveBase.VSFS, 'k:', 'linewidth', 2.5)
            % Quantity cases
            for iQuantity = quantity
                plot(depth, Wave.VSFS{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
            end
            xlabel('Depth (mbsl)')
            ylabel('Shear wave velocity (m/s)')
            axis([axisMinDepth axisMaxDepth 530 600])
            legend( '0 g/dm^3' , '6 g/dm^3' , '15 g/dm^3' , '23 g/dm^3' , '32 g/dm^3' , '40 g/dm^3' )
            title('c')


            %%% Density
            axis3 = subplot(2, 2, 2);
            hold on
            
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, 1.65, 1.8)
            
            % Base case
            plot(depth, WaveBase.bulkDensityFS ./ 1000, 'k:', 'linewidth', 2.5)
            % Quantity cases
            for iQuantity = quantity
                plot(depth, Wave.bulkDensityFS{iQuantity} ./ 1000, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
            end
            xlabel('Depth (mbsl)')
            ylabel('Bulk density (g/cm^3)')
            axis([axisMinDepth axisMaxDepth 1.65 1.8])
            title('b')

            %%% Bulk modulus
            axis4 = subplot(2, 2, 1);
            hold on
            
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, 3e9, 8e9)
            
            % Base case
            plot(depth, WaveBase.bulkKFS, 'k:', 'linewidth', 2.5)
            % Quantity cases
            for iQuantity = quantity
                plot(depth, Wave.bulkKFS{iQuantity}, 'Color', colorStream(iQuantity,:)', 'linewidth', 2.5);
            end
            xlabel('Depth (mbsl)')
            ylabel('Bulk modulus (Pa)')
            axis([axisMinDepth axisMaxDepth 3e9 8e9])
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
            BEQL3P = 481 + obj.seafloorDepth; % mbsf
            plot([BEQL3P, BEQL3P], [-1, 1], '--', 'Color', [.4 .4 .4], 'linewidth', 1.5);
            
            figureCellArray = cell(numel(quantity), 1);
            figureNumber = 0;
            
            for iQuantity = quantity
                figureNumber = figureNumber + 1;
                figureCellArray{figureNumber} = plot(obj.depthArrayFull, Wave.seismogram{iQuantity}, ...
                                                    'Color', colorStream(iQuantity,:)', ...
                                                    'linewidth', 2.5);
            end
%             axis([380 + obj.seafloorDepth, 580 + obj.seafloorDepth, obj.axisMinAmplitude obj.axisMaxAmplitude])
            axis([415 + obj.seafloorDepth, 545 + obj.seafloorDepth, obj.axisMinAmplitude obj.axisMaxAmplitude])
            xlabel('Depth (mbsl)')
            ylabel('Amplitude')
            title('a) Depth Series')
            legend( [figureCellArray{:}], '6 g/dm^3' , '15 g/dm^3' , '23 g/dm^3' , '32 g/dm^3' , '40 g/dm^3' )
            
            
            %%% Time series seismogram
            axis6 = subplot(2, 1, 2);
            hold on
            for iQuantity = quantity
                plot(Wave.time{iQuantity}, Wave.seismogram{iQuantity}, ...
                    'Color', colorStream(iQuantity,:)', ...
                    'Linewidth', 2.5)
            end       
%             axis([4.1 4.3 obj.axisMinAmplitude obj.axisMaxAmplitude])
            axis([4.135 4.28 obj.axisMinAmplitude obj.axisMaxAmplitude])
            xlabel('TWT time (seconds)')
            ylabel('Amplitude')
            title('b) Time Series')
%             legend( '6 g/dm^3' , '15 g/dm^3' , '23 g/dm^3' , '32 g/dm^3' , '40 g/dm^3' )
            
            axis5.Position = [.13 .56 .79 .39];
            axis6.Position = [.13 .07 .79 .39];
            
            figure2.Position(3) = 416;
            set(findall(figure2,'-property','FontSize'),'FontSize',8)
            set(findall(figure2,'-property','FontName'),'FontName','Arial')
        end
        function PlotBackgroundProperties( obj , dataBase )
            
            BSRTop = 456 + obj.seafloorDepth; % mbsl
            BSRBottom = 472 + obj.seafloorDepth; % mbsl
            bulkEQLDepth = 481 + obj.seafloorDepth; % mbsl
            bulkEQLLine = [bulkEQLDepth bulkEQLDepth];
            
            depthAxis = [160 500] + obj.seafloorDepth;
            
            plotLineWidth = 1.5;
            markerLineWidth = 1.25;
            labelLineFont = 8;
            
            figure1 = figure();
            
            figure1.Position(3) = 800;
            set(findall(figure1, '-property', 'FontSize'), 'FontSize', 8)
            set(findall(figure1, '-property', 'FontName'), 'FontName', 'Arial')
            
            interval = 0.13;
            gap = 0.01;
            xStart = 0.06;            
            
            bitRadius = 10.125 ./ 2; % inches
            
            n = 7;
            xCell = cell(n, 1);
            yCell = cell(n, 1);
            limitCell = cell(n, 1);
            axisCell = cell(n, 1);
            titleCell = cell(n, 1);
            xLabelCell = cell(n, 1);
            xTickCell = cell(n, 1);
            
            % Where the tick labels are
            xTickCell{1} = [30 65 100];
            xTickCell{2} = [-10 -5 0 5 10];
            xTickCell{3} = [0.6 1 1.8]; % log scale
            xTickCell{4} = [0.45 0.6 0.75];
            xTickCell{5} = [1.2 1.5 1.8];
            xTickCell{6} = [1500 1700 1900];
            xTickCell{7} = [250 500 750];
            
            % Where the axis limits are
            limitCell{1} = [12 110];
            limitCell{2} = [-13 13];
            limitCell{3} = [0.5 2.1]; % log scale
            limitCell{4} = [0.4 0.8];
            limitCell{5} = [1 2];
            limitCell{6} = [1400 2000];
            limitCell{7} = [100 900];            
            
            xCell{1} = dataBase.GammaRay;
            xCell{2} = obj.CALI.Diameter;
            xCell{3} = dataBase.Resistivity;
            xCell{4} = dataBase.Porosity;
            xCell{5} = dataBase.BulkDensity ./ 1000;
            xCell{6} = dataBase.VP;
            xCell{7} = dataBase.VS;
            
            yCell{1} = dataBase.Depth + obj.seafloorDepth;
            yCell{2} = obj.CALI.Depth + obj.seafloorDepth;
            yCell{3} = dataBase.Depth + obj.seafloorDepth;
            yCell{4} = dataBase.Depth + obj.seafloorDepth;
            yCell{5} = dataBase.Depth + obj.seafloorDepth;
            yCell{6} = dataBase.Depth + obj.seafloorDepth;
            yCell{7} = dataBase.Depth + obj.seafloorDepth;
            
            titleCell{1} = 'a';
            titleCell{2} = 'b';
            titleCell{3} = 'c';
            titleCell{4} = 'd';
            titleCell{5} = 'e';
            titleCell{6} = 'f';
            titleCell{7} = 'g';
            
            xLabelCell{1} = 'Gamma ray (gAPI)';
            xLabelCell{2} = 'Hole radius (in)';
            xLabelCell{3} = 'Resistivity (ohm-m)';
            xLabelCell{4} = 'Porosity';
            xLabelCell{5} = 'Bulk density (g/cm^3)';
            xLabelCell{6} = 'P-velocity (m/s)';
            xLabelCell{7} = 'S-velocity (m/s)';
            
            
            
            for i = 1:n
                
                leftLimit = limitCell{i}(1);
                rightLimit = limitCell{i}(2);    
                xLine = [leftLimit rightLimit];
                
                if i == 7
                    temp = axisCell{6}.Position;
                    axisCell{6}.Position = axisCell{5}.Position;
                end
                
                axisCell{i} = subplot(1, 7, i);
                axisCell{i}.FontSize = 8;
                axisCell{i}.XTick = xTickCell{i};

                
                if i == 7
                    axisCell{6}.Position = temp;
                end                
                
                if i ~= 1
                    axisCell{i}.YTickLabel = [];
                end
                
                if i ~= 3
                    hold on
                end
                
                if i ~= 3
                    DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
                end

                if i == 2
                    caliperRadius = xCell{i} ./ 2;
                    plot(caliperRadius, yCell{i}, 'k', 'Linewidth', 0.5)
                    plot(-caliperRadius, yCell{i}, 'k', 'Linewidth', 0.5)
                    plot([bitRadius bitRadius], depthAxis, 'k-.', 'Linewidth', 1)
                    plot([-bitRadius -bitRadius], depthAxis, 'k-.', 'Linewidth', 1)
                elseif i == 3
                    semilogx(xLine, bulkEQLLine, 'k--', 'Linewidth', markerLineWidth)
                    hold on
                    DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
                    semilogx(xCell{i}, yCell{i}, 'k', 'Linewidth', plotLineWidth)
                    grid on
                    axisCell{i}.XTick = xTickCell{i};
                    axisCell{i}.YTickLabel = [];
                    axisCell{i}.FontSize = 8;
                else
                    plot(xCell{i}, yCell{i}, 'k', 'Linewidth', plotLineWidth)
                end
                
                set(gca, 'YDir', 'Reverse')
                xlabel(xLabelCell{i}, 'Fontsize', 8)
                axis([leftLimit rightLimit depthAxis(1) depthAxis(2)])
                title(titleCell{i})
                
                if i == 1
                    ylabel('Depth (mbsl)')
%                     text(-13, (BSRTop + BSRBottom)/2, 'BSR', 'Fontsize', labelLineFont)
                    text(15, (BSRTop + BSRBottom)/2, 'BSR', 'Fontsize', labelLineFont)
                    text(-28, bulkEQLLine(1), '3P EQL', 'Fontsize', labelLineFont)                    
                end
                
                if i ~= 3
                    plot(xLine, bulkEQLLine, 'k--', 'Linewidth', markerLineWidth)                
                end
                axisCell{i}.Position(1) = xStart;
                axisCell{i}.Position(3) = interval - gap;
                xStart = xStart + interval;
            end
        end
        function PlotThicknessVsQuantity( obj , Wave )
            
            startIndex = 6;
            xData = obj.quantityArray(startIndex : end);
            yData = [Wave.thickness3P{startIndex : end}]';
            
            logFitType = fittype('a + b*log(x)', ...
                'dependent', {'y'}, ...
                'independent', {'x'}, ...
                'coefficients', {'a', 'b'});
            [logFit logGOF] = fit(xData, yData, logFitType)
            
            figure
            axis1 = plot(logFit, 'r-', ...
                        xData, yData, 'ks');

            set(axis1, 'Linewidth', 1.5)
            xlabel('Methane quantity (g/dm^3 of pore volume)')
            ylabel('Transition zone thickness (m)')
            axis([xData(1) xData(end) -inf inf])
            set(gca, 'XScale', 'log')
            legend('Data points', 'Logarithmic fit')
        end
        
        
        function PlotSeismogramOriginalResolution( obj , Wave )
            figure1 = figure();
            startIndex = 6;
            nQuantity = numel(Wave.seismogram);
            
            
            depth = obj.depthArrayFull;
            
            
            %%% Depth series seismogram
            axis1 = subplot(2, 1, 1);
            hold on
            colorStream = colormap(jet(nQuantity));
            
            
            BEQL3P = 481 + obj.seafloorDepth; % mbsl
            plot([BEQL3P, BEQL3P], [-1, 1], '--', 'Color', [.4 .4 .4], 'Linewidth', 1.5)
            
            for iQuantity = startIndex:nQuantity
                plot(depth, Wave.seismogram{iQuantity}, ...
                        'Color', colorStream(iQuantity,:)', ...
                        'Linewidth', 1);
            end
%             axis([380 + obj.seafloorDepth 580 + obj.seafloorDepth -0.15 0.1])
            axis([415 + obj.seafloorDepth, 545 + obj.seafloorDepth, -0.15, 0.1])
            xlabel('Depth (mbsl)')
            ylabel('Amplitude')
            title('a) Depth Series')
            c1 = colorbar('Limits', [6 40], ...
                          'Ticks', linspace(6, 40, 5), ...
                          'TickLabels', linspace(6, 40, 5));
            caxis([6 40])
            c1.Label.String = 'Methane quantity (g/dm^3 of pore volume)';
            c1.Label.FontSize = 8;
            
            %%% Time series seismogram
            axis2 = subplot(2, 1, 2);
            hold on
            colorStream = colormap(jet(nQuantity));
            for iQuantity = startIndex:nQuantity
                plot(Wave.time{iQuantity}, Wave.seismogram{iQuantity}, ...
                    'Color', colorStream(iQuantity,:)', ...
                    'Linewidth', 1)
            end
%             axis([4.15 4.35 -0.15 0.1])
            axis([4.19 4.33 -0.15 0.1])
            xlabel('TWT time (seconds)')
            ylabel('Amplitude')
            title('b) Time Series')
            c2 = colorbar('Limits', [6 40], ...
                          'Ticks', linspace(6, 40, 5), ...
                          'TickLabels', linspace(6, 40, 5));
            caxis([6 40])
            c2.Label.String = 'Methane quantity (g/dm^3 of pore volume)';
            c2.Label.FontSize = 8;
            
            axis1.Position = [.13 .56 .63 .39];
            axis2.Position = [.13 .07 .61 .39];
            
            figure1.Position(3) = 416;
            set(findall(figure1,'-property','FontSize'),'FontSize',8)
            set(findall(figure1,'-property','FontName'),'FontName','Arial')
        end
        function PlotVelocityStructureOriginalResolution( obj , WaveBase , Wave )
            quantity = obj.selectedQuantities;
            depth = obj.depthArray + obj.seafloorDepth;          
            
            figure
            colorStream = jet(numel(Wave.VPFS));
            hold on
            
            leftLimit = 3248;
            rightLimit = 3276;            
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, 1200, 2400)

            
            %%% VP
            % Base case
            plot(depth, WaveBase.VPFS, 'k:', 'Linewidth', 2.5)
            % Methane quantity cases
            for iQuantity = quantity
                
                plot(depth, Wave.VPFS{iQuantity}, ...
                    'Color', colorStream(iQuantity,:)', ...
                    'Linewidth', 2.5);
            end
            xlabel('Depth (mbsl)')
            ylabel('Compressional wave velocity (m/s)')
            axis([450 + obj.seafloorDepth 510 + obj.seafloorDepth 1200 2400])
            legend('0 g/dm^3', '6 g/dm^3', '15 g/dm^3', '23 g/dm^3', '32 g/dm^3', '40 g/dm^3')
            title('Preupscaled Compressional Wave Velocity')
        end
        function PlotPeakAmplitudeRatio( ~ , Wave )
            %%% Smoothing calculation
            % This averages the peak amplitudes based on unique thicknesses 
            % to make it more smooth and understandable
            
            thickness = cell2mat(Wave.thickness3P);
            leadingPeakAmp = Wave.peak(:, 1);
            trailingPeakAmp = Wave.peak(:, 2);
            
            emptyLogical = cellfun(@isempty, leadingPeakAmp);
            
            thickness(emptyLogical) = [];
            leadingPeakAmp(emptyLogical) = [];
            trailingPeakAmp(emptyLogical) = [];
            
            leadingPeakAmp = cell2mat(leadingPeakAmp);
            trailingPeakAmp = cell2mat(trailingPeakAmp);
            
            [uniqueThickness, uniqueIndices, ~] = unique(thickness);
            n = numel(uniqueThickness);
            
            uniqueAvgLeadingPeakAmp = zeros(n, 1);
            uniqueAvgTrailingPeakAmp = zeros(n, 1);
            
            for i = 1:n
                iFrom = uniqueIndices(i);
                if i == n
                    % Average current unique index to end of full data
                    iTo = numel(thickness);
                else
                    % Average current unique index to 1 before next unique
                    % index
                    iTo = uniqueIndices(i + 1) - 1;
                end
                uniqueAvgLeadingPeakAmp(i) = mean(leadingPeakAmp(iFrom : iTo));
                uniqueAvgTrailingPeakAmp(i) = mean(trailingPeakAmp(iFrom : iTo));
            end
            
            uniquePeakAmpRatio = uniqueAvgLeadingPeakAmp ./ uniqueAvgTrailingPeakAmp;
            
            
            figure
            hold on
            
            plot([0 30], [1 1], '--', 'Color', [.4 .4 .4], 'linewidth', 1.5)
            
            plot(uniqueThickness, uniquePeakAmpRatio, 'ks', 'Linewidth', 1.5);
            xlabel('Transition zone thickness (m)')
            ylabel('Leading peak/trailing peak amplitude')
            axis([0 30 -inf inf])
%             axis([0 30 0.8 1.3])
            
            
            %{
            switch caseString
                case 'ParameterSensitivity'
                    exclude = uniqueThickness > 15;
                    axisVector = [0 30 -inf inf];
                case 'OriginalResolution'
                    exclude = uniqueThickness > 15 | uniqueThickness < 3;
                    axisVector = [0 30 -inf inf];
            end
            
            [linearFit, linearGOF] = fit(uniqueThickness, uniquePeakAmpRatio, ...
                                            'poly1', 'Exclude', exclude)
            
            figure
            axis1 = plot(linearFit, 'r-', ...
                        uniqueThickness, uniquePeakAmpRatio, 'ks', ...
                        exclude , 'b*');
            set(axis1, 'LineWidth', 1.5)
            xlabel('Transition zone thickness (m)')
            ylabel('Leading peak/trailing peak amplitude')
            axis(axisVector)
            legend('Data points', 'Excluded data', 'Linear fit')
            %}
        end
        function PlotDickensSeismogram( obj , WaveOriginalResolution , WaveDickens )
            quantity = obj.selectedQuantities;
            colorStream = jet(numel(WaveOriginalResolution.seismogram));
            
            figure1 = figure();
            
            %%% Depth series seismogram
            axis1 = subplot(2, 1, 1);
            hold on
            % 3 phase bulk equilibrium depth line
            BEQL3P = 481 + obj.seafloorDepth; % mbsl
            plot([BEQL3P, BEQL3P], [-1, 1], '--' , 'Color' , [.4 .4 .4] , 'linewidth' , 1.5 );
            
            figureCellArray = cell(numel(quantity) + 1, 1);
            figureNumber = 0;
            
            for iQuantity = quantity
                figureNumber = figureNumber + 1;
                figureCellArray{figureNumber} = plot(obj.depthArrayFull, WaveOriginalResolution.seismogram{iQuantity}, ...
                                                    'Color', colorStream(iQuantity,:)', ...
                                                    'linewidth', 2.5);
            end
            
            figureCellArray{figureNumber + 1} = plot(obj.depthArrayFull, WaveDickens.seismogram, ...
                                                    'k:', ...
                                                    'linewidth', 2.5);
            
            axis([415 + obj.seafloorDepth, 545 + obj.seafloorDepth, -0.15 0.15])
            xlabel('Depth (mbsl)')
            ylabel('Amplitude')
            title('a) Depth Series')
            legend( [figureCellArray{:}], '6 g/dm^3' , '15 g/dm^3' , '23 g/dm^3' , '32 g/dm^3' , '40 g/dm^3' , 'Variable' )
            
            
            %%% Time series seismogram
            axis2 = subplot(2, 1, 2);
            hold on
            for iQuantity = quantity
                plot(WaveOriginalResolution.time{iQuantity}, WaveOriginalResolution.seismogram{iQuantity}, ...
                    'Color', colorStream(iQuantity,:)', ...
                    'Linewidth', 2.5)
            end
            
            plot(WaveDickens.time, WaveDickens.seismogram, ...
                    'k:', ...
                    'Linewidth', 2.5)
            
            axis([4.19 4.33 -0.15 0.15])
            xlabel('TWT time (seconds)')
            ylabel('Amplitude')
            title('b) Time Series')
            
            axis1.Position = [.15 .56 .78 .39];
            axis2.Position = [.15 .07 .78 .39];
            
            figure1.Position(3) = 416;
            figure1.Position(4) = 674; %
            figure1.Position(2) = 10; %
            set(findall(figure1,'-property','FontSize'),'FontSize',8)
            set(findall(figure1,'-property','FontName'),'FontName','Arial')
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
            
            methaneQuantity = Dickens.MethaneQuantity(:, 1);
            depth = Dickens.MethaneQuantity(:,2) + obj.seafloorDepth;
            
            axis_1 = subplot(2,1,1);
            h = plot(methaneQuantity, depth, 'k');
            set(gca,'YDir','Reverse')
            axis([0 40 100 + obj.seafloorDepth 600 + obj.seafloorDepth])
            % grid on
            xlabel('Methane quantity (g/dm^3 of pore volume)')
            ylabel('Depth (mbsl)')
            title('a) Experimentally Measured Methane Quantity')
            set(h,'LineWidth',2.5);
            
            
%             mean(methaneQuantity(depth >= 3248 & depth <= 3272))
            

            % Fluid saturations vs depth
            % figure
            % scatter( Sh_interpolated , Dickens.MethaneQuantity(:,2) )
            % hold on
            % scatter( Sg_interpolated , Dickens.MethaneQuantity(:,2) )
            % hold on
            axis_2 = subplot(2,1,2);
            h = plot( Dickens.Hydrate , obj.depthArray + obj.seafloorDepth );
            set(h,'LineWidth',2.5);
            set(h,'Color',[0 .5 0]);
            hold on
            h = plot( Dickens.Gas , obj.depthArray + obj.seafloorDepth, 'r--' );
            set(h,'LineWidth',2.5);
            set(gca,'YDir','Reverse')
            axis([0 0.2 100 + obj.seafloorDepth 600 + obj.seafloorDepth])
            % axis([0 0.2 420 520])
            xlabel('Fluid saturations')
            ylabel('Depth (mbsl)')
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
        function DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
            rectangle('Position', [leftLimit, BSRTop, rightLimit - leftLimit, BSRBottom - BSRTop], ...
                        'FaceColor', [0.8, 0.8, 0.8], ...
                        'EdgeColor', 'k');            
        end
    end
end























