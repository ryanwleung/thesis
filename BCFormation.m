classdef BCFormation < handle
    properties
        depthArray
        seafloorDepth
        minDepth
        maxDepth
        depthIncrement

        temperatureGradient
        seafloorTemperature
        salinityWtPercent
        methaneMassFractionInHydrate
        
        
        sNwArray
        radiusArray
        lengthArray
        nArray
        
        
        phi0
        phiInf
        B
    end
    properties (Constant)
        waterDensity = 1024;        % kg H2O/m^3 H2O
        hydrateDensity = 928.5;     % kg hydrate/m^3 hydrate
        mwCH4 = 16.043 / 1000;      % kg CH4/mol CH4
        mwH2O = 18.01528 / 1000;    % kg H2O/mol H2O
    end
    methods
        %%% Constructor
        function [ obj ] = BCFormation()
            % Calculations
            obj.methaneMassFractionInHydrate = 4*obj.mwCH4 / (4*obj.mwCH4 + 23*obj.mwH2O);

        end
                
        %%% Main methods
        function [ exportTable , transitionZoneProperties ] = RunSolubilitySaturationRoutine( obj , ch4Quantity )
            
            %%% Get depth and depth dependent parameters
            depth = obj.depthArray;
            pressure = obj.CalcPressure( depth );
            temperature = obj.CalcTemperature( depth );
            gasDensity = BCFormation.CalcGasDensityArray( pressure , temperature );
            
            %%% Load ch4 quantity array
            ch4Quantity = ch4Quantity .* ones(numel(depth), 1);
            % ch4Quantity = ch4Quantity .* zeros(numel(depth), 1);
            % ch4Quantity(1:3:end) = 0;

            %%% Instantiate solubility class object
            solubilityCalculator = BCSolubilityUtil();
            
            %%% Bulk solubility calculations
            [ solBulkLG , solBulkLH , T3PArray ] = ...
                solubilityCalculator.CalcBulkSolubilities( pressure , temperature , obj.salinityWtPercent );
            
            [ ~ , index3PBulkSolEQL ] = min(abs(temperature - T3PArray));
            
            %%% Get bulk saturation from bulk solubilities
            sgBulk = obj.CalcSg2P( ch4Quantity , solBulkLG , gasDensity );
            shBulk = obj.CalcSh2P( ch4Quantity , solBulkLH );

            %%% Get max solubility (smallest pore) for LG and LH
            %%% along with associated 2P saturations
            [ solMaxLG , sg2P ] = obj.CalcMaxSolLG( ch4Quantity , pressure , gasDensity , solBulkLG );
            [ solMaxLH , sh2P ] = obj.CalcMaxSolLH( ch4Quantity , temperature , solBulkLH );

            %%% Get min solubility (largest pore associated with entry
            %%% capillary pressure) for LG and LH
            solMinLG = BCFormation.CalcSolubilityLG( solBulkLG , obj.CalcPcgw(0) * 1e6 , pressure );
            solMinLH = BCFormation.CalcSolubilityLH( solBulkLH , obj.CalcPchw(0) * 1e6 , temperature );

            %%% Get transition zone characteristics (top, bottom, 3P zone
            %%% thickness, and 3P indices)
            %%%
            %%% NOTE: returned bottom3PIndex is estimated using max LG and
            %%% min LH
            top3PIndex = BCFormation.GetTop3PIndex( solMinLG , solMaxLH );
            [ bottom3PIndex , ~ , ~ ] = ...
                BCFormation.GetBottom3PIndex( solMaxLG , solMinLH , obj.depthArray );
            thickness = BCFormation.GetThickness3PZone( depth , top3PIndex , bottom3PIndex );
            indexArrayOf3PZone = top3PIndex : bottom3PIndex;
            
            %%% Setting sg above 3P zone BOTTOM -> 0
            %%%            below 3P zone BOTTOM -> sg 2P
            sg = sg2P;
            sg(1:bottom3PIndex) = 0;
            
            %%% Setting sh above 3P zone TOP -> sh 2P
            %%%            below 3P zone TOP -> 0
            sh = sh2P;
            sh(top3PIndex:end) = 0;
            
            %%% Setting overall solubility to 2P max sol LG and LH cases
            sol = solMaxLH;
            sol(bottom3PIndex + 1:end) = solMaxLG(bottom3PIndex + 1:end);
            
            gasSolubilityPhase2 = solMaxLG;
            
            
            %%% Getting saturations in 3P zone and inserting into arrays
            %%% Returned bottom3PIndex is updated with actual result
            [ sg3P , sh3P , sol3P , bottom3PIndex ] = obj.Calc3P( ch4Quantity , indexArrayOf3PZone , ...
                                                pressure , temperature , gasDensity , ...
                                                solBulkLG , solBulkLH , solMaxLH(top3PIndex) , gasSolubilityPhase2 , sg2P );
            sg(indexArrayOf3PZone) = sg3P;
            sh(indexArrayOf3PZone) = sh3P;
            sol(indexArrayOf3PZone) = sol3P;
            
            
            
            %%% Calculate gas capillary pressures
            %pcgw2PMPa = obj.CalcPcgw(sg2P);
            pcgw2PMPa = obj.CalcPcgw(sgBulk);
            pcgw2PPa = pcgw2PMPa .* 1e6;
            
            pcgw3PMPa = obj.CalcPcgw(sg);
            pcgw3PPa = pcgw3PMPa .* 1e6;
            
            
            
            %%% Compiling arrays into a table for exporting
            exportTable = table();
            exportTable.Depth = depth;
            exportTable.Pressure = pressure;
            exportTable.Temperature = temperature;
            exportTable.Temperature3P = T3PArray;
            exportTable.GasDensity = gasDensity;
            exportTable.GasBulkSol = solBulkLG;
            exportTable.GasMinSol = solMinLG;
            exportTable.GasMaxSol = solMaxLG;
            exportTable.HydrateBulkSol = solBulkLH;
            exportTable.HydrateMinSol = solMinLH;
            exportTable.HydrateMaxSol = solMaxLH;
            %exportTable.GasSat2P = sg2P;
            exportTable.GasSat2P = sgBulk;
            %exportTable.HydrateSat2P = sh2P;
            exportTable.HydrateSat2P = shBulk;
            exportTable.GasSat3P = sg;
            exportTable.HydrateSat3P = sh;
            exportTable.OverallSol = sol;
            exportTable.Pcgw2PPa = pcgw2PPa;
            exportTable.Pcgw3PPa = pcgw3PPa;
            
            %%% Compiling transition zone properties into a structure
            transitionZoneProperties.Top3PIndex = top3PIndex;
            transitionZoneProperties.Bottom3PIndex = bottom3PIndex;
            transitionZoneProperties.Thickness = [];
            transitionZoneProperties.Bulk3PSolEQLIndex = index3PBulkSolEQL;
        end
        function [ exportTable ] = RunRockAndRatioRoutine( obj , exportTable )
            [exportTable.bulkDensityKg, exportTable.porosity] = obj.EstimateBulkDensity();
            exportTable.rockStrengthPa = obj.CalcRockStrength(exportTable);
            
            
            rockStrengthPa = exportTable.rockStrengthPa;
            pcgw2PPa = exportTable.Pcgw2PPa;
            pcgw3PPa = exportTable.Pcgw3PPa;
            
            exportTable.ratio2P = pcgw2PPa ./ rockStrengthPa;
            exportTable.ratio3P = pcgw3PPa ./ rockStrengthPa;            
        end        
        
        
        %%% General calculations
        function [ pressure ] = CalcPressure( obj , depth )
            hydrostaticGradient = obj.waterDensity / 1000 * .433 .* 22620.6; % pa/m (converted from psi/ft)
            pressure = hydrostaticGradient.*(depth + obj.seafloorDepth); % Pa

        end
        function [ temperature ] = CalcTemperature( obj , depth )
            % obj.temperatureGradient is in C deg/km
            temperatureGradientMeters = obj.temperatureGradient ./ 1000; % C deg/m
            temperature = obj.seafloorTemperature + depth .* temperatureGradientMeters; % C
            temperature = temperature + 273.15; % K
        end
        function [ bulkDensityKg , porosity ] = EstimateBulkDensity( obj )
            % Including the non-logged depths (in mbsf) in the effective vertical stress
            % This function is only used for the fracture code
            
            rhoFluid = 1.024; % g/cc, seawater
            rhoGrain = 2.7;   % g/cc, smectite
            
            porosity = obj.phiInf + (obj.phi0 - obj.phiInf) .* exp(-obj.depthArray ./ obj.B);
            bulkDensityG = porosity .* rhoFluid + (1 - porosity) .* rhoGrain;
            bulkDensityKg = bulkDensityG .* 1000;
        end
        function [ rockStrengthPa ] = CalcRockStrength( obj , exportTable )
            % Need to change bulkDensity from g to kg
            verticalEffectiveStress = (exportTable.bulkDensityKg - obj.waterDensity )  .* 9.81; % Pa/m
            minHorizontalEffectiveStress = (obj.poissonsRatio / (1 - obj.poissonsRatio)) .* verticalEffectiveStress;
            
            rockStrengthPa = cumsum(minHorizontalEffectiveStress .* obj.depthIncrement);
            
            
        end

        %%% 2P calculations
        % Saturations
        function [ gasSaturation ] = CalcSg2P( obj , ch4Quantity , solubilityLG , gasDensity )
            % Input: Mass_CH4_Quantity [kg/m^3], Sol_Liq_Gas [mol CH4/kg H2O],
            % Density_Brine [g H2O/cm^3], Density_Gas [kg CH4/m^3]
            % Output: Sg [dimensionless volume fraction of PV]
            kgCH4InCubicMeterWater = solubilityLG .* obj.mwCH4 .* obj.waterDensity;
            
            
            gasSaturation = (ch4Quantity - kgCH4InCubicMeterWater) ./ ...
                            (gasDensity - kgCH4InCubicMeterWater);
                            
            gasSaturation(ch4Quantity < kgCH4InCubicMeterWater) = 0;
        end
        function [ hydrateSaturation ] = CalcSh2P( obj , ch4Quantity , solubilityLH )
            kgCH4InCubicMeterWater = solubilityLH .* obj.mwCH4 .* obj.waterDensity;

            hydrateSaturation = (ch4Quantity - kgCH4InCubicMeterWater) ./ ...
                                (obj.hydrateDensity .* obj.methaneMassFractionInHydrate - kgCH4InCubicMeterWater);
            
            hydrateSaturation(ch4Quantity < kgCH4InCubicMeterWater) = 0;
        end
        % LG solubility
        function [ maxSolLG , sg2P ] = CalcMaxSolLG( obj , ch4Quantity , pressure , gasDensity , gasBulkSolubility )
            n = numel(pressure);
            tempSg2P = zeros(n, 1);
            tempSolLG2P = zeros(n, 1);
            
            gasSaturationBulk2P = obj.CalcSg2P( ch4Quantity , gasBulkSolubility , gasDensity );
            
            for i = 1:n
                %%% If calculated bulk sg is < 0, then all methane is dissolved
                if gasSaturationBulk2P(i) <= 0
                    tempSg2P(i) = 0;
                    tempSolLG2P(i) = gasBulkSolubility(i);
                    continue    
                end


                sg = gasSaturationBulk2P(i);
                
                doWhileFlag = true;
                iteration = 0;
                iterationFactor = 1;
                deltaCellArray = cell(1, 1);
                while doWhileFlag || abs(deltaSg) > 1e-7 
                    doWhileFlag = false;
                    iteration = iteration + 1;
                    
                    
                    [solLGIterated, sgIterated] = obj.CalcMaxSolLGIteration(sg, gasBulkSolubility(i), ch4Quantity(i), pressure(i), gasDensity(i));
                    deltaSg = sg - sgIterated;
                    
                    sg = sgIterated;
                    
                    if isnan(sg)
                        error('NaN found when calculating MaxSolLG, sg2P = %.5f', gasSaturationBulk2P(i))
                    end
                    
                    
                    deltaCellArray{1} = deltaSg;
                end
                
                % %%% Print run status
                % BCFormation.PrintIterationData( 'CalcMaxSolLG' , i , n , iteration , iterationFactor , deltaCellArray )
                
                tempSg2P(i) = sg;
                tempSolLG2P(i) = solLGIterated;
            end
            sg2P = tempSg2P;
            maxSolLG = tempSolLG2P;
        end
        function [ solLG , sg ] = CalcMaxSolLGIteration( obj , sg , bulkSolLG , ch4Quantity , pressure , gasDensity )
            % all scalar input parameters
            
            pcgwMPa = obj.CalcPcgw( sg );
            pcgwPa = pcgwMPa .* 1e6; % convert from MPa to Pa
            
            solLG = BCFormation.CalcSolubilityLG( bulkSolLG , pcgwPa , pressure );
            
            sg = obj.CalcSg2P( ch4Quantity , solLG , gasDensity );
        end
        % LH solubility
        function [ maxSolLH , sh2P ] = CalcMaxSolLH( obj , ch4Quantity , temperature , hydrateBulkSolubility )
            n = numel(temperature);
            tempSh2P = zeros(n, 1);
            tempSolLH2P = zeros(n, 1);
            
            hydrateSaturationBulk2P = obj.CalcSh2P(ch4Quantity, hydrateBulkSolubility);

            for i = 1:n
                if hydrateSaturationBulk2P(i) <= 0
                    tempSh2P(i) = 0;
                    tempSolLH2P(i) = hydrateBulkSolubility(i);
                    continue    
                end


                sh = hydrateSaturationBulk2P(i);
                
                doWhileFlag = true;
                iteration = 0;
                iterationFactor = 1;
                deltaCellArray = cell(1, 1);
                while doWhileFlag || abs( deltaSh ) > 1e-7
                    doWhileFlag = false;
                    iteration = iteration + 1;
                    
                    [solLHIterated, shIterated] = obj.CalcMaxSolLHIteration(sh, hydrateBulkSolubility(i), ch4Quantity(i), temperature(i));
                    deltaSh = sh - shIterated;
                    
                    sh = shIterated;

                    if isnan(sh)
                        error('NaN found when calculating MaxSolLH, sh2P = %.5f', hydrateSaturationBulk2P(i))
                    end
                    
                    
                    deltaCellArray{1} = deltaSh;
                end
                
                % %%% Print run status
                % BCFormation.PrintIterationData( 'CalcMaxSolLH' , i , n , iteration , iterationFactor , deltaCellArray )
                
                tempSh2P(i) = sh;
                tempSolLH2P(i) = solLHIterated;
            end
            
            
            
            sh2P = tempSh2P;
            maxSolLH = tempSolLH2P;
        end        
        function [ solLH , sh ] = CalcMaxSolLHIteration( obj , sh , bulkSolLH , ch4Quantity , temperature )
            % all scalar input parameters
            
            pchwMPa = obj.CalcPchw( sh );
            pchwPa = pchwMPa .* 1e6; % convert from MPa to Pa

            solLH = BCFormation.CalcSolubilityLH(bulkSolLH, pchwPa, temperature);
            sh = obj.CalcSh2P(ch4Quantity, solLH);            
        end
        
        %%% 3P calculations
        function [ sg3P , sh3P , adjustedSol , bottom3PIndex ] = Calc3P( obj , ch4Quantity , indexArrayOf3PZone , ...
                                                            pressure , temperature , gasDensity , ...
                                                            gasBulkSolubility , hydrateBulkSolubility , ...
                                                            hydrateMaxSolubilityAtTop , solubilityPhase2 , gasSaturation2P )
            n = numel(indexArrayOf3PZone);
            sg3P = zeros(n, 1);
            sh3P = zeros(n, 1);
            adjustedSol = zeros(n, 1);
            
            reached2ndPhase = false;
            bottom3PIndex = [];
            %%% Start of Newton's method
            
            % Perturbation for slope calculation
            eps = 1e-5;
            
            % Initial guess of solubility
            solubility = hydrateMaxSolubilityAtTop;
            
            % figure
            % hold on
            
            i = 1;
            while i <= n
                i3P = indexArrayOf3PZone(i);
                
                if reached2ndPhase
                    
                    solubility = solubilityPhase2(i3P);
                    
                    sh = 0;
                    sg = gasSaturation2P(i3P);
                else
                    % Get previous Sg for inital guess of Sg
                    if i == 1
                        sg = 0;
                    else
                        sg = sg3P(i - 1);
                    end
                    % Do while loop for Newton's method
                    % Condition is when the LG and LH solubilities become equal
                    doWhileFlag = true;
                    iteration = 0;
                    iterationFactor = 1;
                    deltaCellArray = cell(1, 1);
                    while doWhileFlag || abs(deltaSol) > 1e-6
                        doWhileFlag = false;
                        iteration = iteration + 1;
                        
                        % Calculating f(x)
                        [ sh , solubilityLG , solubilityLH ] = Calc3PNewtonIteration( obj , ch4Quantity(i3P) , ...
                                                                                    sg , solubility , ...
                                                                                    pressure(i3P) , temperature(i3P) , gasDensity(i3P) , ...
                                                                                    gasBulkSolubility(i3P) , hydrateBulkSolubility(i3P) );
                        deltaSol = solubilityLG - solubilityLH;
                        
                        % Calculating f'(x)
                        [ ~ , solubilityLGPerturbed , solubilityLHPerturbed ] = Calc3PNewtonIteration( obj , ch4Quantity(i3P) , ...
                                                                                    sg + eps , solubility , ...
                                                                                    pressure(i3P) , temperature(i3P) , gasDensity(i3P) , ...
                                                                                    gasBulkSolubility(i3P) , hydrateBulkSolubility(i3P) );
                        deltaSolPerturbed = solubilityLGPerturbed - solubilityLHPerturbed;
                        slope = (deltaSolPerturbed - deltaSol)/eps;
                        
                        % Calculating next iteration sg
                        sg = sg - iterationFactor * deltaSol/slope;
                        % Update solubility by taking the average of the 2
                        solubility = (solubilityLG + solubilityLH)/2;
                        
                        
                        if isnan(sg) || isnan(sh) || isnan(solubilityLG) || isnan(solubilityLH)
                            error('NaN found when calculating 3P saturations')
                        end
                        
                        deltaCellArray{1} = deltaSol;
                        if mod(iteration, 20) == 0 && iterationFactor > 0.01
                            iterationFactor = iterationFactor * 0.9;
                            BCFormation.PrintIterationData( 'Calc3P' , i , n , iteration , iterationFactor , deltaCellArray )
                        end
                    end
                    
                    if solubility > solubilityPhase2(i3P)
                        if reached2ndPhase
                            disp('------------------------ 2nd phase of 3P calc activated twice ------------------------')
                            error('2nd phase of 3P calc activated twice')
                        end
                        reached2ndPhase = true;
                        bottom3PIndex = i3P;
                        continue
                    end
                end                
                % pcgw = BCFormation.CalcPcgwFromSolLG( gasBulkSolubility(i3P) , solubility , pressure(i3P) );
                % pchw = BCFormation.CalcPchwFromSolLH( hydrateBulkSolubility(i3P) , solubility , temperature(i3P) );
                % radiusG = BCFormation.CalcRadiusGasFromPcgw( pcgw );
                % radiusH = BCFormation.CalcRadiusHydrateFromPchw( pchw );
                % scatter(radiusG, obj.depthArray(i3P), 'r', 'filled')
                % scatter(radiusH, obj.depthArray(i3P), 'g', 'filled')
                
                sg3P(i) = sg;
                sh3P(i) = sh;
                adjustedSol(i) = solubility;
                
                i = i + 1;
            end
            % xlabel('Pore size in m^3')
            % ylabel('Depth in mbsf')
            % set(gca, 'YDir', 'Reverse')
        end
        function [ sh , solubilityLG , solubilityLH ] = Calc3PNewtonIteration( obj , ch4Quantity , ...
                                                                                sg , solubility , ...
                                                                                pressure , temperature , gasDensity , ...
                                                                                gasBulkSolubility , hydrateBulkSolubility )
            % Note: sh will go negative if guessed sg is too high to match
            % the input CH4 quantity (since gas = CH4 mass, and CH4 is
            % dissolved in the water
            
            kgCH4InCubicMeterWater = solubility * obj.mwCH4 * obj.waterDensity;
            
            sh = ( ch4Quantity - kgCH4InCubicMeterWater + sg * (kgCH4InCubicMeterWater - gasDensity) ) / ...
                    (obj.hydrateDensity * obj.methaneMassFractionInHydrate - kgCH4InCubicMeterWater );
            
            sg(ch4Quantity < kgCH4InCubicMeterWater) = 0;
            sh(ch4Quantity < kgCH4InCubicMeterWater) = 0;

            test = 2;
            switch test
                case 1
                    % Fraction of phase in a pore vs hydrate
                    competitionFractionOfGas = 0.5;
                    competitionFractionOfHydrate = 1 - competitionFractionOfGas;

                    % Above equal pore size invasion point
                    if( sh/(sh + sg) > competitionFractionOfHydrate)
                        adjustedSg = sg / competitionFractionOfGas;
                        adjustedSh = sh + sg;
                    % Below equal pore size invasion point
                    else
                        adjustedSg = sg + sh;
                        adjustedSh = sh / competitionFractionOfHydrate;
                    end
                case 2
                    % Hydrate in the smallest pores
                    adjustedSg = sg;
                    adjustedSh = sh + sg;
                case 3
                    % Gas in the smallest pores
                    adjustedSg = sg + sh;
                    adjustedSh = sh;
            end
            
            
            pcgwMPa = obj.CalcPcgw( adjustedSg );
            pchwMPa = obj.CalcPchw( adjustedSh );
            
            pcgwPa = pcgwMPa .* 1e6; % convert from MPa to Pa
            pchwPa = pchwMPa .* 1e6; % convert from MPa to Pa
            
            
            solubilityLG = BCFormation.CalcSolubilityLG( gasBulkSolubility , pcgwPa , pressure );
            solubilityLH = BCFormation.CalcSolubilityLH( hydrateBulkSolubility , pchwPa , temperature );
        end
        
        %%% Unused calculations
        function [ sg , sh ] = Calc3P2ndPhase( obj , targetSolubility , shInitialGuess , pressure , temperature , bulkSolubilityLG , bulkSolubilityLH  , ...
                                                    ch4Quantity , gasDensity )
            eps = -1e-6;
            shGuess = shInitialGuess;
            gramsCH4InCubicMeterWater = targetSolubility * obj.mwCH4 * obj.waterDensity;
            
            
            
            pcgw = BCFormation.CalcPcgwFromSolLG( bulkSolubilityLG , targetSolubility , pressure );
            pchw = BCFormation.CalcPchwFromSolLH( bulkSolubilityLH , targetSolubility , temperature );
            
            radiusG = BCFormation.CalcRadiusGasFromPcgw( pcgw );
            radiusH = BCFormation.CalcRadiusHydrateFromPchw( pchw );
            if radiusH < radiusG
                error('radiusH is less than radiusG when it should be larger')
            end
            
            
            
            
            
            % obj.PlotCumPSD();
            % hold on
            
            
            
            
            
            % Do while loop for Newton's method
            % Condition is when the LG and LH solubilities become equal
            doWhileFlag = true;
            iteration = 0;
            iterationFactor = 1;
            while doWhileFlag || abs(deltaSh) > 1e-6
                doWhileFlag = false;
                iteration = iteration + 1;
                
                
                
                
                
                
                
                
                
                % Calculating f(x)
                [ sgFromCumPSD , radiusGPrime ] = obj.InterpCumPSD( shGuess , radiusG , radiusH );
                
                shMassBal = ( ch4Quantity - gramsCH4InCubicMeterWater + sgFromCumPSD * (gramsCH4InCubicMeterWater - gasDensity) ) / ...
                        (obj.hydrateDensity * obj.methaneMassFractionInHydrate - gramsCH4InCubicMeterWater );                
                
                
                deltaSh = shMassBal - shGuess
                
                
                
                % Calculating f'(x)
                [ sgFromCumPSDPerturbed , ~ ] = obj.InterpCumPSD( shGuess + eps , radiusG , radiusH );
                
                shMassBalPerturbed = ( ch4Quantity - gramsCH4InCubicMeterWater + sgFromCumPSDPerturbed * (gramsCH4InCubicMeterWater - gasDensity) ) / ...
                        (obj.hydrateDensity * obj.methaneMassFractionInHydrate - gramsCH4InCubicMeterWater );                
                
                
                deltaShPerturbed = shMassBalPerturbed - (shGuess + eps);
                
                slope =   (deltaShPerturbed - deltaSh) ...
                        / (eps);
                
                
                
                
                % Calculating next iteration shGuess
                shGuess = shGuess - iterationFactor * deltaSh/slope;
                
                
                
                
                if isnan(sgFromCumPSD) || isnan(shGuess)
                    error('NaN found when calculating phase 2 stuff')
                end
            end
            
            sh = shGuess;
            sg = sgFromCumPSD;
            
        end
        
        %%% Plotting methods
        function GenerateResultPlots( obj , exportTable , transitionZoneProperties )
            lineStyle2D = cell(1,3);
            lineStyle2D{1} = 'r--';
            lineStyle2D{2} = 'g--';
            
            lineStyle3D = cell(1,3);
            lineStyle3D{1} = 'r-';
            lineStyle3D{2} = 'g-';
            
            lineStylePc = cell(1,3);
            lineStylePc{1} = 'r--';
            lineStylePc{2} = 'r-';            
            
            lineStyleRatio = cell(1,3);
            lineStyleRatio{1} = 'r--';
            lineStyleRatio{2} = 'r-';    
            
            
            
            solFigure = figure();
            sat2PFigure = figure();
            sat3PFigure = figure();
            pcgwFigure = figure();
            ratioFigure = figure();

            solFigure = obj.PlotSol( solFigure , exportTable , true );
            
            sat2PFigure = obj.PlotSat2P( sat2PFigure , exportTable , transitionZoneProperties , lineStyle2D );
            sat3PFigure = obj.PlotSat3P( sat3PFigure , exportTable , lineStyle3D );
            
            pcgwFigure = obj.PlotRockStrength( pcgwFigure  , exportTable );
            pcgwFigure = obj.PlotPcgw( pcgwFigure , exportTable , transitionZoneProperties , lineStylePc );
            ratioFigure = obj.PlotRatio( ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio );
            
        end
        function [ solFigure ] = PlotSol( ~ , solFigure , exportTable , doPlotBulkAndMinSol )
            depth = exportTable.Depth;
            solBulkLG = exportTable.GasBulkSol;
            solMinLG = exportTable.GasMinSol;
            solMaxLG = exportTable.GasMaxSol;
            solBulkLH = exportTable.HydrateBulkSol;
            solMinLH = exportTable.HydrateMinSol;
            solMaxLH = exportTable.HydrateMaxSol;
            sol = exportTable.OverallSol;
            
            width = 2;
            
            figure(solFigure)
            hold on
            if doPlotBulkAndMinSol
                plot( solBulkLG , depth , 'r-' , 'linewidth' , width )
                plot( solBulkLH , depth , 'g-' , 'linewidth' , width )
                plot( solMinLG , depth , 'r-.' , 'linewidth' , width )
                plot( solMinLH , depth , 'g-.' , 'linewidth' , width )
            end
            plot( solMaxLG , depth , 'r--' , 'linewidth' , width )
            plot( solMaxLH , depth , 'g--' , 'linewidth' , width )
            plot( sol , depth , 'b' , 'linewidth' , width )
            
            xlabel('CH_4 Solubility (mol CH4/kg H2O)')
            ylabel('Depth (mbsf)')
            set(gca, 'YDir', 'Reverse')
            legend('Bulk G-W', 'Bulk H-W', 'Min G-W', 'Min H-W', 'Max G-W', 'Max H-W', 'Calculated solubility')
        end       
        function [ sat2PFigure ] = PlotSat2P( ~ , sat2PFigure , exportTable , transitionZoneProperties , lineStyle )
            bulkEquilibrium3PIndex = transitionZoneProperties.Bulk3PSolEQLIndex;
            depth = exportTable.Depth;
            
            sg2P = exportTable.GasSat2P;
            sg2P(1 : bulkEquilibrium3PIndex - 1) = 0;
            
            sh2P = exportTable.HydrateSat2P;
            sh2P(bulkEquilibrium3PIndex : end) = 0;     
            
            [ depthSg , sg2P ] = BCFormation.GetModifiedPlotArrays2P( depth , sg2P , bulkEquilibrium3PIndex );
            [ depthSh , sh2P ] = BCFormation.GetModifiedPlotArrays2P( depth , sh2P , bulkEquilibrium3PIndex );

            figure(sat2PFigure)
            hold on
            plot( sg2P , depthSg , lineStyle{1} , 'linewidth' , 3 )
            plot( sh2P , depthSh , lineStyle{2} , 'linewidth' , 3 )
            xlabel('Saturation')
            ylabel('Depth (mbsf)')
            legend('Gas', 'Hydrate')
            set(gca, 'YDir', 'Reverse')
        end
        function [ sg3PFigure ] = PlotSat3P( ~ , sg3PFigure , exportTable , lineStyle )
            depth = exportTable.Depth;
            sg3P = exportTable.GasSat3P;
            sh3P = exportTable.HydrateSat3P;
            
            figure(sg3PFigure)
            hold on
            plot( sg3P , depth , lineStyle{1} , 'linewidth' , 3 )
            plot( sh3P , depth , lineStyle{2} , 'linewidth' , 3 )
            xlabel('Saturation')
            ylabel('Depth (mbsf)')
            legend('Gas', 'Hydrate')
            set(gca, 'YDir', 'Reverse')
        end
        function [ pcgwFigure ] = PlotRockStrength( obj , pcgwFigure , exportTable )
            
            depth = obj.depthArray;
            rockStrengthPa = exportTable.rockStrengthPa;
            
            figure(pcgwFigure)
            
            plot(rockStrengthPa ./ 1e6, depth, 'k', 'linewidth', 3)
            hold on
            
            xlabel('Pressure (MPa)')
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
        end
        function [ pcgwFigure ] = PlotPcgw( obj , pcgwFigure , exportTable , transitionZoneProperties , lineStylePc )
            
            depth = obj.depthArray;
            pcgw2PPa = exportTable.Pcgw2PPa;
            pcgw3PPa = exportTable.Pcgw3PPa;
            
            %%% Hotfix to eleiminate gas capillary entry pressure
            zeroSgLogical = exportTable.GasSat3P == 0;
            pcgw3PPa(zeroSgLogical) = 0;
            indexToSmooth = find(~zeroSgLogical, 1);
            [ depthFor3P , pcgw3PPa ] = BCFormation.GetModifiedPlotArrays2P( depth , pcgw3PPa , indexToSmooth );
            
            
            bulkEquilibrium3PIndex = transitionZoneProperties.Bulk3PSolEQLIndex;
            pcgw2PPa(1 : bulkEquilibrium3PIndex - 1) = 0;
            [ depthFor2P , pcgw2PPa ] = BCFormation.GetModifiedPlotArrays2P( depth , pcgw2PPa , bulkEquilibrium3PIndex );
            
            
            
            
            figure(pcgwFigure)
            hold on
            
            plot( pcgw2PPa ./ 1e6 , depthFor2P , lineStylePc{1} , 'linewidth' , 3 )
            plot( pcgw3PPa ./ 1e6 , depthFor3P , lineStylePc{2} , 'linewidth' , 3 )
            
            xlabel('Pressure (MPa)')
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
            legend('Minimum horizontal effective stress', 'Bulk equilibrium model', 'Three-phase stability model')
        end
        function [ ratioFigure ] = PlotRatio( obj , ratioFigure , exportTable , transitionZoneProperties , lineStyleRatio )
            
            depth = obj.depthArray;
            rockStrengthPa = exportTable.rockStrengthPa;
            pcgw2PPa = exportTable.Pcgw2PPa;
            pcgw3PPa = exportTable.Pcgw3PPa;
            
                        
            ratio2P = pcgw2PPa ./ rockStrengthPa;
            ratio3P = pcgw3PPa ./ rockStrengthPa;
            
            bulkEquilibrium3PIndex = transitionZoneProperties.Bulk3PSolEQLIndex;
            ratio2P(1 : bulkEquilibrium3PIndex - 1) = 0;
            [ depthFor2P , ratio2P ] = BCFormation.GetModifiedPlotArrays2P( depth , ratio2P , bulkEquilibrium3PIndex );
            
            %%% Hotfix to eleiminate gas capillary entry pressure
            zeroSgLogical = exportTable.GasSat3P == 0;
            ratio3P(zeroSgLogical) = 0;
            indexToSmooth = find(~zeroSgLogical, 1);
            [ depthFor3P , ratio3P ] = BCFormation.GetModifiedPlotArrays2P( depth , ratio3P , indexToSmooth );

            
            figure(ratioFigure)
            hold on
            
            plot( [1 1] , [obj.minDepth obj.maxDepth] , 'k--' , 'linewidth' , 2 , 'HandleVisibility','off' )
            
            plot( ratio2P , depthFor2P , lineStyleRatio{1} , 'linewidth' , 3 )
            plot( ratio3P , depthFor3P , lineStyleRatio{2} , 'linewidth' , 3 )            
            
            xlabel('Gas overpressure/minimum horizontal effective stress')
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
            legend('Bulk equilibrium model', 'Three-phase stability model')
        end
    end
    methods (Static)
        %%% Utility
        function [ outputArray ] = ArrayInsert( baseArray , topDepthIndex , inputArray )
            
            outputArray = baseArray;
            bottomDepthIndex = topDepthIndex + length(inputArray) - 1;
            
            outputArray( topDepthIndex : bottomDepthIndex ) = inputArray;
            
        end
        function [ depth , variable ] = GetModifiedPlotArrays2P( depth , variable , iInsert )
            
            depth = [ depth(1 : iInsert - 1) ; depth(iInsert) ; depth(iInsert : end) ];
            variable = [ variable(1 : iInsert - 1) ; variable(iInsert - 1) ; variable(iInsert : end) ];
            
        end
        
        function PrintIterationData( functionString , i , n , iteration , iterationFactor , deltaCellArray )
            
            if mod(i, 25) ~= 0 && i ~= 1 && i ~= n
                return
            end
            
            disp('--------------------------------------------------------')
            disp(['Solving: ', functionString])
            disp(['Step ', num2str(i), ' of ', num2str(n), ': (', num2str(round(i / n * 100)), '%)']) 
            disp(['Iteration: ', num2str(iteration)])
            disp(['Iteration factor: ', num2str(iterationFactor)])
            nDelta = numel(deltaCellArray);
            for iDelta = 1:nDelta
                disp(['Delta ', num2str(iDelta), ': ', num2str(deltaCellArray{iDelta})])
            end
        end
        
        %%% Get 3P zone properties
        function [ top3PIndex ] = GetTop3PIndex( gasMinSolubility , hydrateMaxSolubility )
            difference = hydrateMaxSolubility - gasMinSolubility;
            difference(difference < 0) = inf;
            [ ~ , top3PIndex ] = min(difference);
        end
        function [ bottom3PIndex , gasMaxSolubilityAtBottom , actualBottom3PDepth ] = GetBottom3PIndex( gasMaxSolubility , hydrateMinSolubility , depthArray )
            difference = gasMaxSolubility - hydrateMinSolubility;
            difference(difference < 0) = inf;
            [ ~ , bottom3PIndex ] = min(difference);
            
            % gasMaxSolubility line
            sol1 = [gasMaxSolubility(bottom3PIndex) gasMaxSolubility(bottom3PIndex + 1)];
            depth1 = [depthArray(bottom3PIndex) depthArray(bottom3PIndex + 1)];
            
            % hydrateMinSolubility line
            sol2 = [hydrateMinSolubility(bottom3PIndex) hydrateMinSolubility(bottom3PIndex + 1)];
            depth2 = [depthArray(bottom3PIndex) depthArray(bottom3PIndex + 1)];
            
            % polynomial order 1 for intersection
            p1 = polyfit(sol1, depth1, 1);
            p2 = polyfit(sol2, depth2, 1);
            
            % find intersection
            gasMaxSolubilityAtBottom = fzero(@(x) polyval(p1 - p2, x), 3);
            actualBottom3PDepth = polyval(p1, gasMaxSolubilityAtBottom);
        end
        function [ thickness ] = GetThickness3PZone( depth , top3PIndex , bottom3PIndex )
            thickness = abs( depth(top3PIndex) - depth(bottom3PIndex) );
        end
        
        %%% Get gas density using PR1978
        function [ gasDensityArray] = CalcGasDensityArray( pressureArray , temperatureArray )
            nCount = length(pressureArray);
            gasDensityArray = zeros(nCount,1);
            
            pr1978Calculator = BCPR1978Util();
            for iCount = 1:nCount
                [ gasDensity ] = pr1978Calculator.CalcGasDensity( pressureArray(iCount) , temperatureArray(iCount) );
                gasDensityArray(iCount) = gasDensity;
                
            end
        end

        %%% Calculate LG and LH solubility by adjusting bulk solubility
        %%% with capillary pressure
        function [ solubilityLG ] = CalcSolubilityLG( bulkSolubilityLG , pcgw , waterPressure )
            solubilityLG = bulkSolubilityLG + bulkSolubilityLG .* ( pcgw ./ waterPressure );
        end
        function [ solubilityLH ] = CalcSolubilityLH( bulkSolubilityLH , pchw , temperature )
            % bulkSolLH units don't matter (this is a ratio calculation)
            % standard solubility units are mol CH4 / kg H2O
            % pchw in Pa
            % temperature in K
            
            hydrateStoichiometryFactor = 5.75;
            hydrateLatticeWaterMolarVolume = 22.6; % cm^3 / mol            
            hydrateLatticeWaterMolarVolume = hydrateLatticeWaterMolarVolume ./ 100^3; % cm^3/mol -> m^3/mol
            RGasConstant = 8.3144598; % J / mol K = Pa m^3 / mol K

            solubilityLH = bulkSolubilityLH + bulkSolubilityLH .* ( pchw .* hydrateStoichiometryFactor .* hydrateLatticeWaterMolarVolume ...
                                                                    ./ (RGasConstant .* temperature) );
        end
        
        
        
        function [ pcgw ] = CalcPcgwFromSolLG( bulkSolubilityLG , solubilityLG , waterPressure )
            pcgw = waterPressure .* (solubilityLG - bulkSolubilityLG) ./ bulkSolubilityLG;
        end
        function [ pchw ] = CalcPchwFromSolLH( bulkSolubilityLH , solubilityLH , temperature )
            % bulkSolLH units don't matter (this is a ratio calculation)
            % standard solubility units are mol CH4 / kg H2O
            % pchw in Pa
            % temperature in K
            
            hydrateStoichiometryFactor = 5.75;
            hydrateLatticeWaterMolarVolume = 22.6; % cm^3 / mol            
            hydrateLatticeWaterMolarVolume = hydrateLatticeWaterMolarVolume ./ 100^3; % cm^3/mol -> m^3/mol
            RGasConstant = 8.3144598; % J / mol K = Pa m^3 / mol K

            pchw = (solubilityLH - bulkSolubilityLH) ./ bulkSolubilityLH ...
                    ./ ( hydrateStoichiometryFactor .* hydrateLatticeWaterMolarVolume ...
                        ./ (RGasConstant .* temperature) );
        end
        
        function [ radiusG ] = CalcRadiusGasFromPcgw( pcgw )
            % pcgw in Pa
            % radiusG in m^2
            
            radiusG = 2 * 0.072 ./ pcgw;
        end
        function [ radiusH ] = CalcRadiusHydrateFromPchw( pchw )
            % pchw in Pa
            % radiusH in m^2
            
            radiusH = 2 * 0.027 ./ pchw;
        end
    end

end


























