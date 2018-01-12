classdef BCFormation < handle
    properties
        Bulk % remove this later

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
        
        
        
    end
    properties (Constant)
        waterDensity = 1024;        % kg H2O/m^3 H2O
        hydrateDensity = 900;       % kg hydrate/m^3 hydrate
        poissonsRatio = 0.4;    
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

            %%% Instantiate solubility class object
            solubilityCalculator = BCSolubilityUtil();
            
            %%% Bulk solubility calculations
            [ solBulkLG , solBulkLH , T3PArray ] = ...
                solubilityCalculator.CalcBulkSolubilities( pressure , temperature , obj.salinityWtPercent );
            
            [ ~ , index3PBulkSolEQL ] = min( abs(temperature - T3PArray) );
            
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
            top3PIndex = BCFormation.GetTop3PIndex( solMinLG , solMaxLH );
            bottom3PIndex = BCFormation.GetBottom3PIndex( solMaxLG , solMinLH );
            thickness = BCFormation.GetThickness3PZone( depth , top3PIndex , bottom3PIndex );
            indexArrayOf3PZone = top3PIndex + 1:bottom3PIndex - 1;
            
            %%% Setting sg above 3P zone BOTTOM -> 0
            %%%            below 3P zone BOTTOM -> sg 2P
            sg = sg2P;
            sg(1:bottom3PIndex - 1) = 0;
            
            %%% Setting sh above 3P zone TOP -> sh 2P
            %%%            below 3P zone TOP -> 0
            sh = sh2P;
            sh(top3PIndex + 1:end) = 0;

            %%% Setting overall solubility to 2P max sol LG and LH cases
            sol = solMaxLH;
            sol(bottom3PIndex:end) = solMaxLG(bottom3PIndex:end);
            
            %%% Getting saturations in 3P zone and inserting into arrays
            [ sg3P , sh3P , sol3P ] = obj.Calc3P( ch4Quantity , indexArrayOf3PZone , ...
                                                pressure , temperature , gasDensity , ...
                                                solBulkLG , solBulkLH , solMaxLH(top3PIndex) );
            sg(indexArrayOf3PZone) = sg3P;
            sh(indexArrayOf3PZone) = sh3P;
            sol(indexArrayOf3PZone) = sol3P;

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
            exportTable.GasSat2P = sg2P;
            exportTable.HydrateSat2P = sh2P;
            exportTable.GasSat3P = sg;
            exportTable.HydrateSat3P = sh;
            exportTable.OverallSol = sol;
            
            %%% Compiling transition zone properties into a structure
            transitionZoneProperties.Top3PIndex = top3PIndex;
            transitionZoneProperties.Bottom3PIndex = bottom3PIndex;
            transitionZoneProperties.Thickness = thickness;
            transitionZoneProperties.Bulk3PSolEQLIndex = index3PBulkSolEQL;
        end
        
        %%% General calculations
        function [ pressure ] = CalcPressure( obj , depth )
            hydrostaticGradient = obj.waterDensity / 1000 * .433 .* 22620.6; % pa/m (converted from psi/ft)
%             hydrostaticGradient = 10 * 1e6 / 1000;
            pressure = hydrostaticGradient.*(depth + obj.seafloorDepth); % Pa

        end
        function [ temperature ] = CalcTemperature( obj , depth )
            % obj.temperatureGradient is in C deg/km
            temperatureGradientMeters = obj.temperatureGradient ./ 1000; % C deg/m
            temperature = obj.seafloorTemperature + depth .* temperatureGradientMeters; % C
            temperature = temperature + 273.15; % K
        end
        function [ rockStrength ] = CalcRockStrength( obj )
            verticalEffectiveStress = (obj.DataTable.bulkDensity - obj.waterDensity) .* 1000 .* 9.81 ./ 1e6;
            minHorizontalEffectiveStress = (obj.poissonsRatio / (1 - obj.poissonsRatio)) .* verticalEffectiveStress;
            rockStrength = cumsum(minHorizontalEffectiveStress);
        end        
        
        %%% 2P calculations
        % Saturations
        function [ gasSaturation ] = CalcSg2P( obj , ch4Quantity , solubilityLG , gasDensity )
            % Input: Mass_CH4_Quantity [g/decimeter^3], Sol_Liq_Gas [mol CH4/kg H2O],
            % Density_Brine [g H2O/cm^3], Density_Gas [kg CH4/m^3]
            % Output: Sg [dimensionless volume fraction of PV]
            gramsCH4InCubicMeterWater = solubilityLG .* obj.mwCH4 .* obj.waterDensity;
            
            gasSaturation = (ch4Quantity - gramsCH4InCubicMeterWater) ./ ...
                            (gasDensity - gramsCH4InCubicMeterWater);
            
%             gasSaturation( isnan(gasSaturation) ) = 0; 
            
        end
        function [ hydrateSaturation ] = CalcSh2P( obj , ch4Quantity , solubilityLH )
            gramsCH4InCubicMeterWater = solubilityLH .* obj.mwCH4 .* obj.waterDensity;

            hydrateSaturation = (ch4Quantity - gramsCH4InCubicMeterWater) ./ ...
                                (obj.hydrateDensity .* obj.methaneMassFractionInHydrate - gramsCH4InCubicMeterWater);
            
%             hydrateSaturation( isnan(hydrateSaturation) ) = 0; 
        end
        % LG solubility
        function [ maxSolLG , sg2P ] = CalcMaxSolLG( obj , ch4Quantity , pressure , gasDensity , gasBulkSolubility )
            n = numel(pressure);
            tempSg2P = zeros(n, 1);
            tempSolLG2P = zeros(n, 1);
            
            gasSaturationBulk2P = obj.CalcSg2P( ch4Quantity , gasBulkSolubility , gasDensity );
            
            for i = 1:n
                
                sg = gasSaturationBulk2P(i);
                
                doWhileFlag = true;
                iteration = 0;
                iterationFactor = 0.75;
                deltaCellArray = cell(1, 1);
                while doWhileFlag || abs(deltaSg) > 1e-7 
                    doWhileFlag = false;
                    iteration = iteration + 1;
                    
                    
                    [solLGIterated, sgIterated] = obj.CalcMaxSolLGIteration(sg, gasBulkSolubility(i), ch4Quantity, pressure(i), gasDensity(i));
                    deltaSg = sg - sgIterated;
                    
                    
                    sg = sg - iterationFactor * (sg - sgIterated);
                    
                    deltaCellArray{1} = deltaSg;
                end
                
                BCFormation.PrintIterationData( 'CalcMaxSolLG' , i , n , iteration , iterationFactor , deltaCellArray )
                
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
            
            hydrateSaturationBulk2P = obj.CalcSh2P(ch4Quantity ,hydrateBulkSolubility);

            for i = 1:n
                sh = hydrateSaturationBulk2P(i);
                
                doWhileFlag = true;
                iteration = 0;
                iterationFactor = 0.75;
                deltaCellArray = cell(1, 1);
                while doWhileFlag || abs( deltaSh ) > 1e-7
                    doWhileFlag = false;
                    iteration = iteration + 1;
                    
                    [solLHIterated, shIterated] = obj.CalcMaxSolLHIteration(sh, hydrateBulkSolubility(i), ch4Quantity, temperature(i));
                    deltaSh = sh - shIterated;
                    
                    sh = sh - iterationFactor * (sh - shIterated);
                    
                    deltaCellArray{1} = deltaSh;
                end
                
                BCFormation.PrintIterationData( 'CalcMaxSolLH' , i , n , iteration , iterationFactor , deltaCellArray )
                
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
        function [ sg3P , sh3P , adjustedSol ] = Calc3P( obj , ch4Quantity , indexArrayOf3PZone , ...
                                                        pressure , temperature , gasDensity , ...
                                                        gasBulkSolubility , hydrateBulkSolubility , hydrateMaxSolubilityAtTop )
            n = numel(indexArrayOf3PZone);
            sg3P = zeros(n, 1);
            sh3P = zeros(n, 1);
            adjustedSol = zeros(n, 1);
            
            
            %%% Start of Newton's method
            
            % Perturbation for slope calculation
            eps = 1e-5;
            
            % Initial guess of solubility
            solubility = hydrateMaxSolubilityAtTop;
            
            for i = 1:n
                i3P = indexArrayOf3PZone(i);
                
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
                    [ sh , solubilityLG , solubilityLH ] = Calc3PNewtonIteration( obj , ch4Quantity , ...
                                                                                sg , solubility , ...
                                                                                pressure(i3P) , temperature(i3P) , gasDensity(i3P) , ...
                                                                                gasBulkSolubility(i3P) , hydrateBulkSolubility(i3P) );
                    deltaSol = solubilityLG - solubilityLH;
                    
                    
                    % Calculating f'(x)
                    [ ~ , solubilityLGPerturbed , solubilityLHPerturbed ] = Calc3PNewtonIteration( obj , ch4Quantity , ...
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
                
                
                
                sg3P(i) = sg;
                sh3P(i) = sh;
                adjustedSol(i) = solubility;
            end
        end
        function [ sh , solubilityLG , solubilityLH ] = Calc3PNewtonIteration( obj , ch4Quantity , ...
                                                                                sg , solubility , ...
                                                                                pressure , temperature , gasDensity , ...
                                                                                gasBulkSolubility , hydrateBulkSolubility )
            % Note: sh will go negative if guessed sg is too high to match
            % the input CH4 quantity (since gas = CH4 mass, and CH4 is
            % dissolved in the water
            
            gramsCH4InCubicMeterWater = solubility * obj.mwCH4 * obj.waterDensity;
            
            sh = ( ch4Quantity - gramsCH4InCubicMeterWater + sg * (gramsCH4InCubicMeterWater - gasDensity) ) / ...
                    (obj.hydrateDensity * obj.methaneMassFractionInHydrate - gramsCH4InCubicMeterWater );
            
            % Fraction of phase in a pore vs hydrate
%             competitionFractionOfGas = 0.019;
%             competitionFractionOfGas = 0.2;
%             competitionFractionOfGas = 0.35;
%             competitionFractionOfGas = 0.65;
%             competitionFractionOfGas = 0.8;
%             competitionFractionOfGas = 0.981;
            
            competitionFractionOfGas = 0.5;
            competitionFractionOfHydrate = 1 - competitionFractionOfGas;
            
%             % Above equal pore size invasion point
%             if( sh/(sh + sg) > competitionFractionOfHydrate)
%                 adjustedSg = sg / competitionFractionOfGas;
%                 adjustedSh = sh + sg;
%             % Below equal pore size invasion point
%             else
%                 adjustedSg = sg + sh;
%                 adjustedSh = sh / competitionFractionOfHydrate;
%             end
            
            adjustedSg = sg;
            adjustedSh = sh + sg;
            
            pcgwMPa = obj.CalcPcgw( adjustedSg );
            pchwMPa = obj.CalcPchw( adjustedSh );
            
            pcgwPa = pcgwMPa .* 1e6; % convert from MPa to Pa
            pchwPa = pchwMPa .* 1e6; % convert from MPa to Pa
            
            
            solubilityLG = BCFormation.CalcSolubilityLG( gasBulkSolubility , pcgwPa , pressure );
            solubilityLH = BCFormation.CalcSolubilityLH( hydrateBulkSolubility , pchwPa , temperature );
        end
        
        
        
        
        
        %%% Plotting methods
        function GenerateResultPlots( obj , exportTable , transitionZoneProperties )
            lineStyle2D = cell(1,3);
            lineStyle2D{1} = 'r--';
            lineStyle2D{2} = 'g--';
%             lineStyle2D{3} = 'r--';
            
            lineStyle3D = cell(1,3);
            lineStyle3D{1} = 'r-';
            lineStyle3D{2} = 'g-';
%             lineStyle3D{3} = 'r-';
            
            solFigure = figure();
            sat2PFigure = figure();
            sat3PFigure = figure();
%             pcgwFigure = figure();
%             ratioFigure = figure();

            solFigure = obj.PlotSol( solFigure , exportTable );
            sat2PFigure = obj.PlotSat2P( sat2PFigure , exportTable , transitionZoneProperties , lineStyle2D );
            sat3PFigure = obj.PlotSat3P( sat3PFigure , exportTable , lineStyle3D );
            
%             [ pcgwFigure ] = PlotRockStrength( obj , pcgwFigure );
%             [ ratioFigure ] = PlotFractureRatio( obj , ratioFigure );
%             
%             
%             for iStorage = ch4QuantityToPlot
%                 iLineStyle = iLineStyle + 1;
%                 
%                 [ pcgwFigure ] = PlotPcgw( obj , pcgwFigure , iStorage , lineStyle2D{iLineStyle} , lineStyle3D{iLineStyle} );
%                 [ ratioFigure ] = PlotRatio( obj , ratioFigure , iStorage , lineStyle2D{iLineStyle} , lineStyle3D{iLineStyle} );
%                 
%             end
        end
        function [ solFigure ] = PlotSol( ~ , solFigure , exportTable )
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
            
            plot( solBulkLG , depth , 'r-' , 'linewidth' , width )
            hold on
            plot( solBulkLH , depth , 'g-' , 'linewidth' , width )
            hold on
            plot( solMinLG , depth , 'r-.' , 'linewidth' , width )
            hold on
            plot( solMinLH , depth , 'g-.' , 'linewidth' , width )
            hold on
            plot( solMaxLG , depth , 'r--' , 'linewidth' , width )
            hold on
            plot( solMaxLH , depth , 'g--' , 'linewidth' , width )
            hold on            
            plot( sol , depth , 'b' , 'linewidth' , width )
            
            xlabel('CH4 Solubility (mol CH4/kg H2O')
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
            legend('Bulk LG', 'Bulk LH', 'Min LG', 'Min LH', 'Max LG', 'Max LH', 'Actual solubility')
        end       
        function [ sat2PFigure ] = PlotSat2P( ~ , sat2PFigure , exportTable , transitionZoneProperties , lineStyle )
            bulkEquilibrium3PIndex = transitionZoneProperties.Bulk3PSolEQLIndex;
            depth = exportTable.Depth;
            
            sg2P = exportTable.GasSat2P;
            sg2P(1 : bulkEquilibrium3PIndex - 1) = 0;
            
            sh2P = exportTable.HydrateSat2P;
            sh2P(bulkEquilibrium3PIndex + 1 : end) = 0;     
            
            [ depth , sg2P ] = BCFormation.GetModifiedPlotArrays2P( depth , sg2P , bulkEquilibrium3PIndex );
            [ ~ , sh2P ] = BCFormation.GetModifiedPlotArrays2P( depth , sh2P , bulkEquilibrium3PIndex );

            figure(sat2PFigure)
            plot( sg2P , depth , lineStyle{1} , 'linewidth' , 3 )
            hold on
            plot( sh2P , depth , lineStyle{2} , 'linewidth' , 3 )
            xlabel('Saturation')
            ylabel('Depth (mbsf)')
            legend('Gas', 'Hydrate')
            set(gca,'YDir','Reverse')
        end
        function [ sg3PFigure ] = PlotSat3P( ~ , sg3PFigure , exportTable , lineStyle )
            depth = exportTable.Depth;
            sg3P = exportTable.GasSat3P;
            sh3P = exportTable.HydrateSat3P;
            
            figure(sg3PFigure)
            plot( sg3P , depth , lineStyle{1} , 'linewidth' , 3 )
            hold on
            plot( sh3P , depth , lineStyle{2} , 'linewidth' , 3 )
            xlabel('Saturation')
            ylabel('Depth (mbsf)')
            legend('Gas', 'Hydrate')
            set( gca , 'YDir' , 'Reverse' )
        end
        
        % not done yet below

        function [ pcgwFigure ] = PlotRockStrength( obj , pcgwFigure )
            
            depth = obj.DataTable.depth;
            rockStrength = obj.DataTable.rockStrength;
            
            figure(pcgwFigure)
            
            plot( rockStrength , depth , 'k' , 'linewidth' , 3 )
            hold on

            xlabel('Pressure (MPa)')
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')

        end
        function [ pcgwFigure ] = PlotPcgw( obj , pcgwFigure , iStorage , lineStyle2D , lineStyle3D )
            
            depth = obj.DataTable.depth;
            pcgw2P = obj.Storage.pcgw2P(:,iStorage);
            
            sg2P = obj.Storage.sg2P(:,iStorage);
            bulkEquilibrium3PIndex = find( sg2P , 1 );
            
            [ depthFor2P , pcgw2P ] = BCFormation.GetModifiedPlotArrays2P( depth , pcgw2P , bulkEquilibrium3PIndex );
            
            pcgw3P = obj.Storage.pcgw3P(:,iStorage);
            
            
            figure(pcgwFigure)
            
            plot( pcgw2P , depthFor2P , lineStyle2D , 'linewidth' , 3 )
            hold on
            plot( pcgw3P , depth , lineStyle3D , 'linewidth' , 3 )
            hold on
            
            xlabel('Pressure (MPa)')
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
            
        end
        function [ ratioFigure ] = PlotRatio( obj , ratioFigure , iStorage , lineStyle2D , lineStyle3D )
            
            depth = obj.DataTable.depth;
            ratio2P = obj.Storage.ratio2P(:,iStorage);
            ratio3P = obj.Storage.ratio3P(:,iStorage);
            
            figure(ratioFigure)
            
            plot( ratio2P , depth , lineStyle2D , 'linewidth' , 3 )
            hold on
            plot( ratio3P , depth , lineStyle3D , 'linewidth' , 3 )
            hold on
            
            
            xlabel('Gas overpressure / rock strength')
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
            
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
            [ ~ , top3PIndex ] = min(abs( gasMinSolubility - hydrateMaxSolubility ));
        end
        function [ bottom3PIndex ] = GetBottom3PIndex( gasMaxSolubility , hydrateMinSolubility )
            [ ~ , bottom3PIndex ] = min(abs( gasMaxSolubility - hydrateMinSolubility ));
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
    end
    % UNUSED CODE
    %{
            OLD CODE THAT INTERPOLATES FROM LIU AND FLEMINGS TO COMPARE
        
            obj.DataTable.gasBulkSolubilityInterpolated = interp1( obj.Bulk.Depth , obj.Bulk.Solubility , depth );
        
            % % BR adjustment to match Liu bulk sol results
            % gasBulkSolubility = gasBulkSolubility - .002;
            % hydrateBulkSolubility = hydrateBulkSolubility + .0024;        
        
            interpolatedBulkSolLogical = ~isnan(obj.DataTable.gasBulkSolubilityInterpolated);
            tempInterpBulkSol = obj.DataTable.gasBulkSolubilityInterpolated(interpolatedBulkSolLogical);
            tempInterpDepth = obj.DataTable.depth(interpolatedBulkSolLogical);        
        
            switch class(obj)
                case 'DCHydrateRidge'
                    tempInterpDepth = tempInterpDepth + 790;
                    depth = depth + 790;
                case 'DCBlakeRidge'
                    tempInterpDepth = tempInterpDepth + 2780;
                    depth = depth + 2780;
            end     
   
            plot( tempInterpBulkSol , tempInterpDepth , 'k' , 'LineWidth' , 3 )
            hold on        
        
        
        
            switch class(obj)
                case 'DCHydrateRidge'
                    load('Hydrate Ridge Data\Solubility Plots\BulkSolHR.mat')
                    plot( BulkSolHR(:,1) , BulkSolHR(:,2) )    
                case 'DCBlakeRidge'
                    load('Blake Ridge Data\Solubility Plots\BulkSolBR.mat')
                    plot( BulkSolBR(:,1) , BulkSolBR(:,2) )
            end
        
        
        
            % gasMaxSolubility = gasMaxSolubility + 0.0062;
        
        
        
        
        %}
end


























