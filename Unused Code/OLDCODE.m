
% from BCFormation.m
    % UNUSED CODE
    %{
        function [ maxSolLG , sg2P ] = CalcMaxSolLG( obj , ch4Quantity , pressure , gasDensity , gasBulkSolubility )
            n = numel(pressure);
            tempSg2P = zeros(n, 1);
            tempSolLG2P = zeros(n, 1);
            
            eps = 1e-3;
            
            gasSaturationBulk2P = obj.CalcSg2P( ch4Quantity , gasBulkSolubility , gasDensity );
            
            for i = 1:n
                
                sg = gasSaturationBulk2P(i);
                
                doWhileFlag = true;
                iteration = 0;
                iterationFactor = 1;
                deltaCellArray = cell(1, 1);
                while doWhileFlag || abs(deltaSg) > 1e-7 
                    doWhileFlag = false;
                    iteration = iteration + 1;
                    
                    %%% f(x)
                    [solLGIterated, sgIterated] = obj.CalcMaxSolLGIteration(sg, gasBulkSolubility(i), ch4Quantity, pressure(i), gasDensity(i));
                    deltaSg = sg - sgIterated;
                    
                    
                    %%% f'(x)
%                     [solLGPerturbed, sgPerturbed] = obj.CalcMaxSolLGIteration(sg + eps, gasBulkSolubility(i), ch4Quantity, pressure(i), gasDensity(i));
%                     deltaSgPerturbed = sg - sgPerturbed;
%                     slope = (deltaSgPerturbed - deltaSg)/eps;
%                     
%                     % Calculating next iteration sg
%                     sg = sg - iterationFactor * deltaSg/slope;
                    sg = sgIterated;
                    
                    if isnan(sg)
                        error('NaN found when calculating MaxSolLG, sg2P = %.5f', gasSaturationBulk2P(i))
                    end
                    
                    
                    deltaCellArray{1} = deltaSg;
                end
                
                %%% Print run status
                BCFormation.PrintIterationData( 'CalcMaxSolLG' , i , n , iteration , iterationFactor , deltaCellArray )
                
                tempSg2P(i) = sg;
                tempSolLG2P(i) = solLGIterated;
            end
            sg2P = tempSg2P;
            maxSolLG = tempSolLG2P;
        end
        function [ sg3P , sh3P , adjustedSol ] = Calc3PReverse( obj , ch4Quantity , indexArrayOf3PZone , ...
                                                        pressure , temperature , gasDensity , ...
                                                        gasBulkSolubility , hydrateBulkSolubility , hydrateMaxSolubilityAtTop )
            n = numel(indexArrayOf3PZone);
            sg3P = zeros(n, 1);
            sh3P = zeros(n, 1);
            adjustedSol = zeros(n, 1);
            
            
            %%% Start of Newton's method
            
            % Perturbation for slope calculation
            eps = -1e-5;
            
            % Initial guess of solubility
            solubility = hydrateMaxSolubilityAtTop;
            
            for i = n:-1:1
                i3P = indexArrayOf3PZone(i);
                
                % Get previous Sg for inital guess of Sg
                if i == n
                    sg = 0.6;
%                     sg = 0.12;
                else
                    sg = sg3P(i + 1);
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





%% code from plotting log data for background properties figure

%{

                        
            %%% Gamma ray
            leftLimit = 12;
            rightLimit = 110;    
            xLine = [leftLimit rightLimit];
            
            axis1 = subplot(1,7,1);
            axis1.FontSize = 8;
            axis1.XTick = [30 65 100];            
            hold on
            
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
            
            text(-13, (BSRTop + BSRBottom)/2, 'BSR', 'Fontsize', labelLineFont)
            text(-28, bulkEQLLine(1), '3P EQL', 'Fontsize', labelLineFont)
            
            plot(dataBase.GammaRay, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
            xlabel('Gamma ray (gAPI)','FontSize', 8)
            axis([leftLimit rightLimit 160 500])
            title('a')
            
            plot( xLine , bulkEQLLine , 'k--' , 'linewidth' , markerLineWidth )
            
            
            
            
            %%% Caliper
            leftLimit = -11;
            rightLimit = 11;
            xLine = [leftLimit rightLimit];
            
            caliperRadius = obj.CALI.Diameter ./ 2; % inches
            
            axis2 = subplot(1,7,2);
            axis2.FontSize = 8;
            axis2.YTickLabel = [];
            axis2.XTick = [-10 -5 0 5 10];
            hold on
            
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
            
            plot(caliperRadius, obj.CALI.Depth , 'k' , 'linewidth' , .5 )
            plot(-caliperRadius, obj.CALI.Depth , 'k' , 'linewidth' , .5 )
            plot( [bitRadius bitRadius] , [160 500] , 'k-.' , 'linewidth' , 1 )
            plot( [-bitRadius -bitRadius] , [160 500] , 'k-.' , 'linewidth' , 1 )
            
            set(gca,'YDir','Reverse')
            xlabel('Hole radius (in)')
            axis([leftLimit rightLimit 160 500])
            title('b')    
            
            plot( xLine , bulkEQLLine , 'k--' , 'linewidth' , markerLineWidth )    






            %%% Resisitivity
            leftLimit = 0.5;
            rightLimit = 1.5;
            xLine = [leftLimit rightLimit];

            axis3 = subplot(1,7,3);
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
            hold on

            plot(dataBase.Resistivity, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
            set(gca,'YDir','Reverse')
            xlabel('Resistivity (ohm-m)')
            axis([leftLimit rightLimit 160 500])
            title('c')
            hold on
            plot( xLine , bulkEQLLine , 'k--' , 'linewidth' , markerLineWidth )

            axis3.FontSize = 8;
            axis3.YTickLabel = [];
            axis3.XTick = [0.6 1 1.4];






            %%% Calculated porosity
            leftLimit = 0.4;
            rightLimit = 0.8;
            xLine = [leftLimit rightLimit];

            axis4 = subplot(1,7,4);
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
            hold on    
            plot(dataBase.Porosity, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
            plot( Data.log(:,3) , Data.log(:,1) , 'k' , 'linewidth' , plotLineWidth )
            set(gca,'YDir','Reverse')
            xlabel('Porosity')
            axis([leftLimit rightLimit 160 500])
            title('d')  
            hold on
            plot( xLine , bulkEQLLine , 'k--' , 'linewidth' , markerLineWidth )
            
            axis4.FontSize = 8;
            axis4.YTickLabel = [];
            axis4.XTick = [0.45 0.6 0.75];
            
            
            
            
            % Bulk density
            leftLimit = 1.3;
            rightLimit = 1.9;    
            xLine = [leftLimit rightLimit];

            axis5 = subplot(1,7,5);

            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
            hold on

            xText = leftLimit + 0.02*(rightLimit - leftLimit);
            text(xText, (BSRTop + BSRBottom)/2, 'BSR')
            text(xText, bulkEQLLine(1) + 0.6*(BSRBottom - BSRTop), '3P EQL')    


            plot(dataBase.BulkDensity ./ 1000, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
            ylabel('Depth (mbsf)')
            set(gca,'YDir','Reverse')
            xlabel('Bulk density (g/cm^3)')
            axis([leftLimit rightLimit 160 500])
            title('e')
            hold on
            plot( xLine , bulkEQLLine , 'k--' , 'linewidth' , markerLineWidth )    

            axis5.FontSize = 8;


            % VP
            leftLimit = 1400;
            rightLimit = 2000;    
            xLine = [leftLimit rightLimit];

            axis6 = subplot(1,7,6);
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
            rectangle('Position', [leftLimit, BSRTop, rightLimit - leftLimit, BSRBottom - BSRTop], ...
                        'FaceColor', [0.8, 0.8, 0.8], ...
                        'EdgeColor', 'k');
            hold on    
            plot(dataBase.VP, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
            set(gca,'YDir','Reverse')
            xlabel('P-velocity (m/s)')
            axis([leftLimit rightLimit 160 500])
            title('f')
            hold on
            plot( xLine , bulkEQLLine , 'k--' , 'linewidth' , markerLineWidth )    

            axis6.FontSize = 8;
            axis6.YTickLabel = [];
            axis6.XTick = [1500 1700 1900];



            % VS
            leftLimit = 100;
            rightLimit = 900;    
            xLine = [leftLimit rightLimit];

            axis7 = subplot(1,7,7);
            DCSeismicAnalysisBR.DrawRectangle(leftLimit, rightLimit, BSRTop, BSRBottom)
            hold on
            plot(dataBase.VS, dataBase.Depth, 'k', 'linewidth', plotLineWidth)
            set(gca,'YDir','Reverse')
            xlabel('S-velocity (m/s)')
            axis([leftLimit rightLimit 160 500])
            title('g')
            hold on
            plot( xLine , bulkEQLLine , 'k--' , 'linewidth' , markerLineWidth )

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
%}

%% bulk solubility code from Formation class

%{
        function [ gasBulkSolubility , hydrateBulkSolubility ] = calcBulkSolubilities( obj , pressure , temperature )
            n = numel(pressure);
            
            gasBulkSolubility = zeros(n,1);
            hydrateBulkSolubility = zeros(n,1);
            
            
            gramsNaClper1kgWater = (obj.salinityWt/100)*1000 / (1 - (obj.salinityWt/100));
            molality = gramsNaClper1kgWater / 58.44; % moles NaCl / 1 kg water
%             molality = 0.6;


            for i = 1:n
                [ gasBulkSolubility(i) ] = Formation.calcBulkSolubilityLG( pressure(i) , temperature(i) , molality );
                [ hydrateBulkSolubility(i) ] = Formation.calcBulkSolubilityLH( pressure(i) , temperature(i) , molality );
            end


            
            
            
        end
        %}
        
%{
        % Liquid-Gas Solubility EOS
        function [ bulkSolubilityLG ] = calcBulkSolubilityLG( pressure , temperature , molalityNaCl )
            % Input is pressure (Pa), temperature (K), and NaCl (mol CH4/L H2O) at depth Z
            
            x_CH4 = 1; % mole composition in the vapor phase
            
            Tc = 190.6; % K
            Tr = temperature ./ Tc;
            
            Pc = 46.41; % bars
            pressure = pressure .* 1e-5; % changes P from Pa to bars
            Pr = pressure ./ Pc;
              
            % corrected parameters
            % Interaction parameters for dimensionless chemical potential (mu/RT) CH4
            mu_c = [4.30210345e1;  ... 
                -6.83277221e-2; ...
                -5.68718730e3; ...
                3.56636281e-5;  ... 
                -5.79133791e1; ...
                6.11616662e-3;  ...
                -7.85528103e-4; ...
                -9.42540759e-2; ...
                1.92132040e-2;  ...
                -9.17186899e-06];
            
            
            % Interaction parameters for lambda CH4,Na
            lambda_c = [9.92230792e-2;  ...
                2.57906811e-5;  ...
                0;              ...
                0;              ...
                0;              ...
                0;              ...
                0;              ...
                1.83451402e-2;  ...
                0;              ...
                -8.07196716e-6];
            
            % Interaction parameters for xi CH4,Na,Cl
            xi_c = [-6.23943799e-3; ...
                0;              ...
                0;              ...
                0;              ...
                0;              ...
                0;              ...
                0;              ...
                0;              ...
                0;              ...
                0];
            
            
            

            
            [ mu_par ]      = Formation.calcPitzerParameter( mu_c , pressure , temperature );
            [ lambda_par ]  = Formation.calcPitzerParameter( lambda_c , pressure , temperature );
            [ xi_par ]      = Formation.calcPitzerParameter( xi_c , pressure , temperature );
            
            
%             asdf = CSolubilityUtil();
%             [ muPar , lambdaPar , xiPar ] = asdf.CalcPitzerParameters( pressure , temperature );
%             if(muPar - mu_par ~= 0)
%                 disp('mu parameterization not correct')
%                 asdfmu - mu_par
%             end
%             if(lambdaPar - lambda_par ~= 0)
%                 disp('lambda parameterization not correct')
%                 asdflambda - lambda_par
%             end
%             if(xiPar - xi_par ~= 0)
%                 disp('xi parameterization not correct')
%                 asdfxi - xi_par
%             end
            
            
            
            
            
            
            [ Vr , B , C , D , E , G ] = Formation.calcVr(Pr,Tr);
            Z = Pr*Vr/Tr;
            
            fugacityCoefficient = exp( Z - 1 - log(Z) + B/Vr + C/(2*Vr^2) + D/(4*Vr^4) + E/(5*Vr^5) + G );
            
            bulkSolubilityLG = x_CH4 .* pressure / exp( mu_par + 2 * lambda_par * molalityNaCl + xi_par * molalityNaCl^2 - log(fugacityCoefficient) );

        end
        function [ par ] = calcPitzerParameter( c , P , T )
            % Input 1-D array with all 10 interaction parameters in order
            % Output resulting summation of parameters
            % Vectorization already supported
            
            % EQ 7 in Duan et al. (1992)
            par =   c(1) ...
                + c(2).*T ...
                + c(3)./T ...
                + c(4).*T.^2 ...
                + c(5)./(680 - T) ...
                + c(6).*P ...
                + c(7).*P.*log(T) ...
                + c(8).*P./T ...
                + c(9).*P/(680 - T) ...
                + c(10).*(P.^2)./T;
            
        end
        function [ Vr , B , C , D , E , G ] = calcVr( Pr , Tr )
            % Newton's Method iteration to calculate real Vr using definition Vr = (V*Pc)/(R*Tc)
            % Only take scalar Pr and Tr, no vectorization
            
            eps = 1e-4; % perturbation for slope calculation
            
            % Initial guess of Vr using IG assumption for Z = 1
            Vr = Tr/Pr;
            [ deltaZ ] = Formation.calcDuanEOS( Vr , Pr , Tr );
            iteration = 0;
            while( abs(deltaZ) > 1e-10 )
                iteration = iteration + 1;
                
                if(iteration > 20)
                    disp('Iteration exceeded limit of 20 times for Vr calculation for gas bulk solubility calculation')
                    break
                end
                
                dDeltaZ_dVr = (Formation.calcDuanEOS( Vr + eps , Pr , Tr ) - deltaZ) ./ eps;
                
                Vr_new = Vr - (deltaZ ./ dDeltaZ_dVr);
                
                Vr = Vr_new;
                
                [ deltaZ , B , C , D , E , G ] = Formation.calcDuanEOS( Vr , Pr , Tr );
                
            end
        end
        function [ deltaZ , B , C , D , E , G ] = calcDuanEOS( Vr , Pr , Tr )
            % Calculates the delta Z (compressibility factor) by the 2 definitions of
            % Z in Duan et al. (1992)'s appendix. The correct Z and Vr is found when
            % delta Z is returned as ~0
            
            a1 =    8.72553928e-2;
            a2 =    -7.52599476e-1;
            a3 =    3.75419887e-1;
            a4 =    1.07291342e-2;
            a5 =    5.49626360e-3;
            a6 =	-1.84772802e-2;
            a7 =	3.18993183e-4;
            a8 =	2.11079375e-4;
            a9 =	2.01682801e-5;
            a10 =	-1.65606189e-5;
            a11 =	1.19614546e-4;
            a12 =	-1.08087289e-4;
            alpha =	4.48262295e-2;
            beta =	7.53970000e-1;
            gamma =	7.7167000e-2;
            
            B = a1 + a2/(Tr^2) + a3/(Tr^3);
            C = a4 + a5/(Tr^2) + a6/(Tr^3);
            D = a7 + a8/(Tr^2) + a9/(Tr^3);
            E = a10 + a11/(Tr^2) + a12/(Tr^3);
            F = alpha/(Tr^3);
            G = F/(2*gamma)*(beta + 1 - (beta + 1 + gamma/(Vr^2))*exp(-gamma/(Vr^2)));
            
            deltaZ = (Pr*Vr/Tr) - ( 1 + B/Vr + C/(Vr^2) + D/(Vr^4) + E/(Vr^5) + F/(Vr^2)*(beta + gamma/(Vr^2))*exp(-gamma/(Vr^2)) );
        end
        
        
        % Liquid-Hydrate Solubility EOS
        function [ bulkSolubilityLH ] = calcBulkSolubilityLH( pressure , temperature , molalityNaCl ) 
%             R=83.1441;
            
            c1 = 4.30210345e1;
            c2 = -6.83277221e-2;
            c3 = -5.68718730e3;
            c4 = 3.56636281e-5;
            c5 = -5.79133791e1;
            c6 = 6.11616662e-3;
            c7 = -7.85528103e-4;
            c8 = -9.42540759e-2;
            c9 = 1.92132040e-2;
            c10 = -9.17186899e-06;
   
            lambda_c1 = 9.92230792e-2;
            lambda_c2 = 2.57906811e-5;
            lambda_c3 = 0.0;
            lambda_c4 = 0.0;
            lambda_c5 = 0.0;
            lambda_c6 = 0.0;
            lambda_c7 = 0.0;
            lambda_c8 = 1.83451402e-2;
            lambda_c9 = 0.0;
            lambda_c10 = -8.07196716e-6;
            
            xi_c1 = -6.23943799e-3;
            xi_c2 = 0.0;
            xi_c3 = 0.0;
            xi_c4 = 0.0;
            xi_c5 = 0.0;
            xi_c6 = 0.0;
            xi_c7 = 0.0;
            xi_c8 = 0.0;
            xi_c9 = 0.0;
            xi_c10 = 0.0;
            

            P = pressure .* 1e-5; % P converted from Pa to bars 

            T = temperature;
            
            
            lambda_par = lambda_c1 + lambda_c2*T + lambda_c3/T + lambda_c4*T^2 + lambda_c5/(680-T) + lambda_c6*P + lambda_c7*P*log(T) + lambda_c8*P/T + lambda_c9*P/(680-T) + lambda_c10*P^2/T;
            
            xi_par = xi_c1 + xi_c2*T + xi_c3/T + xi_c4*T^2 + xi_c5/(680-T) + xi_c6*P + xi_c7*P*log(T) + xi_c8*P/T + xi_c9*P/(680-T) + xi_c10*P^2/T;
            
            
            mu_par = c1 + c2*T + c3/T + c4*T^2 + c5/(680-T) + c6*P + c7*P*log(T) + c8*P/T + c9*P/(680-T) + c10*P^2/T;
            
            asdf = CSolubilityUtil();
            [ muPar , lambdaPar , xiPar ] = asdf.CalcPitzerParameters( P , T );
            if(muPar - mu_par ~= 0)
                disp('mu parameterization not correct')
                asdfmu - mu_par
            end
            if(lambdaPar - lambda_par ~= 0)
                disp('lambda parameterization not correct')
                asdflambda - lambda_par
            end
            if(xiPar - xi_par ~= 0)
                disp('xi parameterization not correct')
                asdfxi - xi_par
            end
            
            waterActivity = Formation.calcWaterActivity( molalityNaCl );
            
            %fugacity=fsolve('fugacity_methane',0.5*P,optimset('Display','on'),P0,T0,aw);
            
            fugacity = Formation.calcFugacity( pressure , temperature , waterActivity );
            
            bulkSolubilityLH = fugacity*exp(-mu_par - 2 .* lambda_par .* molalityNaCl - xi_par .* molalityNaCl.^2);

        end
        function [ waterActivity ] = calcWaterActivity( molalityNaCl )
            
            A_phi = 0.3915; % at temperature of 25 deg C
            b = 1.2; % unit is kg1/2 / mol1/2
            
            m = molalityNaCl;
            
            mu_M = 1;
            mu_X = 1;
            mu = 2;
            
            z_M = +1;
            z_X = -1;
            
            I = 1/2 .* (m*z_M^2 + m*z_X^2);
            %I=0.6196;
            
            f_phi = -A_phi .* sqrt(I) ./ (1 + b .* sqrt(I));
            
            alpha = 2.0; % unit is kg1/2 / mol1/2
            B0 = 0.0765;
            B1 = 0.2664;
            
            B_MX = B0 + B1 .* exp(-alpha .* sqrt(I));
            
            C_MX = 0.00127;
            
            phi = 1 + abs(z_M*z_X) .* f_phi + m .* (2 .* mu_M .* mu_X ./ mu) .* B_MX + m.^2 .* (2 .* (sqrt(mu_M .* mu_X)).^3 ./ mu) .* C_MX;
            
            waterActivity = exp(-18 .* mu .* m ./ 1000 .* phi);
        end
        function [ fugacityX ] = calcFugacity( pressure , temperature , waterActivity )
            
            fugacityX = 0.5 .* pressure .* 1e-5; % P converted from Pa to bars 
            eps = 1e-4; % perturbation for slope calculation

            [ fugacityY ] = Formation.calcFugacityEOS( fugacityX , pressure , temperature , waterActivity );
            iteration = 0;
            while( abs( fugacityY ) > 1e-10 )
                iteration = iteration + 1;
                
                if(iteration > 20)
                    disp('Iteration exceeded limit of 20 times for fugacity calculation for hydrate bulk solubility calculation')
                    break
                end
                
                
                perturbedFugacityY = Formation.calcFugacityEOS( fugacityX + eps , pressure , temperature , waterActivity );
                dy_dx = (perturbedFugacityY - fugacityY) ./ eps;
                
                
                fugacityXNew = fugacityX - (fugacityY ./ dy_dx);
                
                fugacityX = fugacityXNew;
                
                [ fugacityY ] = Formation.calcFugacityEOS( fugacityXNew , pressure , temperature , waterActivity );
                
            end
        end
        function [ fugacity ] = calcFugacityEOS( x , pressure , temperature , waterActivity )
            
            R = 83.1441; % in ??? / mol K (Gas constant)
            
            
            
            
            % vi is the number of cavitis of type i per molecular of water
            % 46 water, 2 small cavity and 6 large cavity per unit cell.
            v1 = 2/46; % Munck et al 1988 - Table 1
            v2 = 6/46; % Munck et al 1988 - Table 1
            
            % Small cavity
            A1 = 0.7228e-3; % K/atm % Munck et al 1988 - Table 6
            B1 = 3187; % K % Munck et al 1988 - Table 6
            
            % large cavity
            A2 = 23.35e-3; % K/atm % Munck et al 1988 - Table 6
            B2 = 2653; % K % Munck et al 1988 - Table 6
            
            
            delta_mu0 = 12640; % Munck et al 1988 - Table 3
            delta_H0 = -48580; % Munck et al 1988 - Table 3
            delta_V0 = 4.6; % Munck et al 1988 - Table 3
            delta_Cp = 391.6; % Munck et al 1988 - Table 3
            
            P = pressure .* 1e-5;
%             P = pressure;
            T = temperature;
            T_bar = (T + 273.15)/2; % Munck et al 1988 - Equation 8
            
            
            V_beta = 22.6;
            V_l = 18;
            
            r = 20E-9;
            gamma_hw = 0.027;
            
            
            
            left = delta_mu0 ./ (R*273.15) ... % Munck et al 1988 - Equation 7, 8 &&&& Henry 1999 - Equation A4 (water activity term)
                - (delta_H0 - delta_Cp .* 273.15) ./ R .* (-1/T + 1/273.15) - delta_Cp ./ R .* (log(T) - log(273.15)) ...
                + delta_V0 ./ (R .* T_bar) .* P ...
                - log(waterActivity);
            
            
            % + ...
            %V_beta*2*gamma_hw/r/1e5/(R*T);
            
            % left=delta_mu0/(R*273.15) - ...
            %      delta_H0/R*(-1/T+1/273.15) + ...
            %      delta_V0/(R*T)*P;
            
            C1 = A1 ./ T .* exp(B1 ./ T); % Munck et al 1988 - Equation 4
            C2 = A2 ./ T .* exp(B2 ./ T); % Munck et al 1988 - Equation 4
            
            fugacity = left + v1 .* log(1 - C1 .* x ./ (1 + C1 .* x)) + v2 .* log(1 - C2 .* x ./ (1 + C2 .* x));
        end
        %}


















































