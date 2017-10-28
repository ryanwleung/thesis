















% bulk solubility code from Formation class

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


















































