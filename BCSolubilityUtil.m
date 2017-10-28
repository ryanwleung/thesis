classdef BCSolubilityUtil < handle
    properties
        molalityNaCl
        x_NaCl
        
        pitzerConstantsMu_ChemPotential_CH4Liquid
        pitzerConstantsLambda_2ndOrderInteraction_CH4Na
        pitzerConstantsXi_3rdOrderInteraction_CH4NaCl
        
        duanLKConstants
        
        shibueCriticalConstants
        shibueGConstants
        shibueHConstants
        shibueWaterPhiConstants
        
        wagnerPrussDensityConstants
        
        langmuirParameters
        
        cavityParameters
        
        daigleParameters
    end
    properties (Constant)
        % pitzerModel can be:
        % 1992 using Duan (1992) OR
        % 2006 using Duan (2006)
        pitzerModel = 1992;
        
        davieCorrectionFactor = 15

        gasConstantR = 83.14472;    % bars cm^3 / mol K
        
        mwNaCl = 58.44;         % g NaCl/mol CH4
        mwH2O = 18.01528;       % g H2O/mol H2O
        
        
        waterTc = 647.096;      % Critical temperature of pure water in K
                                % See Shibue (2003) for usage in equation 8
        waterPc = 22.064;       % Critical pressure of pure water in MPa
        
        wagnerPrussReferenceDensity = 322;  % kg H2O/m^3 H2O
        
        partialMolarVolumeAqueousCH4 = 35; % cm^3/mole
        
        methaneTc = 190.6;      % K
        methanePc = 46.41;      % bars
        
    end
    methods
        %%% Constructor and setup
        function [ obj ] = BCSolubilityUtil()
            obj.SetPitzerParConstants();
            obj.SetDuanLKConstants();
            obj.SetShibueConstants();
            obj.SetWagnerPrussConstants();
            obj.SetLangmuirParameters();
            obj.SetCavityParameters();
            obj.SetDaigleParameters();
        end
        function SetPitzerParConstants( obj )
            switch obj.pitzerModel
                case 1992
                    % For Duan (1992) "The prediction of..." equation 7 from Table 2
                    obj.pitzerConstantsMu_ChemPotential_CH4Liquid = ...
                        [    4.30210345e1;  ... 4.30210345e1 correct 
                            -6.83277221e-2; ...
                            -5.68718730e3;  ...
                             3.56636281e-5; ... 3.56636281e-5 correct
                            -5.79133791e1;  ...
                             6.11616662e-3; ...
                            -7.85528103e-4; ...
                            -9.42540759e-2; ...
                             1.92132040e-2; ...
                            -9.17186899e-06 ];            
                    obj.pitzerConstantsLambda_2ndOrderInteraction_CH4Na = ...
                        [    9.92230792e-2; ...
                             2.57906811e-5; ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             1.83451402e-2; ...
                             0;             ...
                            -8.07196716e-6 ];  
                    obj.pitzerConstantsXi_3rdOrderInteraction_CH4NaCl = ...
                        [   -6.23943799e-3; ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0 ];
                case 2006
                    % For Duan and Mao (2006) equation 9 from Table 3
                    obj.pitzerConstantsMu_ChemPotential_CH4Liquid = ...
                        [    0.83143711e1;  ... 
                            -0.72772168e-3; ...
                             0.21489858e4;  ...
                            -0.14019672e-4; ... 
                            -0.66743449e6;  ...
                             0.76985890e-2; ...
                            -0.50253331e-5; ...
                            -0.30092013e1;  ...
                             0.48468502e3;  ...
                             0 ];            
                    obj.pitzerConstantsLambda_2ndOrderInteraction_CH4Na = ...
                        [   -0.81222036;    ...
                             0.10635172e-2; ...
                             0.18894036e3;  ...
                             0;             ...
                             0;             ...
                             0.44105635e-4; ...
                             0;             ...
                             0;             ...
                             0;             ...
                            -0.46797718e-10 ];  
                    obj.pitzerConstantsXi_3rdOrderInteraction_CH4NaCl = ...
                        [   -0.29903571e-2; ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0;             ...
                             0 ];
            end
        end
        function SetDuanLKConstants( obj )
            % For Duan and Mao (2006) equation A1 & A4 
            % From Duan (1992) "The prediction of..." Table A1
            obj.duanLKConstants = zeros(15,1);
            obj.duanLKConstants(1)  =    8.72553928e-2;
            obj.duanLKConstants(2)  =   -7.52599476e-1;
            obj.duanLKConstants(3)  =    3.75419887e-1;
            obj.duanLKConstants(4)  =    1.07291342e-2;
            obj.duanLKConstants(5)  =    5.49626360e-3;
            obj.duanLKConstants(6)  =	-1.84772802e-2;
            obj.duanLKConstants(7)  =	 3.18993183e-4;
            obj.duanLKConstants(8)  =	 2.11079375e-4;
            obj.duanLKConstants(9)  =	 2.01682801e-5;
            obj.duanLKConstants(10) =	-1.65606189e-5;
            obj.duanLKConstants(11) =	 1.19614546e-4;
            obj.duanLKConstants(12) =	-1.08087289e-4;
            obj.duanLKConstants(13) =	 4.48262295e-2;
            obj.duanLKConstants(14) =	 7.53970000e-1;
            obj.duanLKConstants(15) =	 7.71670000e-2;
        end
        function SetShibueConstants( obj )
            % For Shibue (2003) euqation 3 and 4
            obj.shibueCriticalConstants = zeros(10,1);
            obj.shibueCriticalConstants(1)  =  8.78054e1;
            obj.shibueCriticalConstants(2)  =  2.42541e3;
            obj.shibueCriticalConstants(3)  = -6.07779e3;
            obj.shibueCriticalConstants(4)  =  1.17033e6;
            obj.shibueCriticalConstants(5)  =  9.00404e2;
            obj.shibueCriticalConstants(6)  = -2.92542e4;
            obj.shibueCriticalConstants(7)  =  1.39806e6;
            obj.shibueCriticalConstants(8)  = -2.80756e7;
            obj.shibueCriticalConstants(9)  =  2.41637e8;
            obj.shibueCriticalConstants(10) = -7.18726e8;
            
            % For Shibue (2003) equation 8
            obj.shibueGConstants = zeros(6,1);
            obj.shibueGConstants(1) = -7.85951783;
            obj.shibueGConstants(2) =  1.84408259;
            obj.shibueGConstants(3) = -11.7866497;
            obj.shibueGConstants(4) =  22.6807411;
            obj.shibueGConstants(5) = -15.9618719;
            obj.shibueGConstants(6) =  1.80122502;
            
            % For Shibue (2003) equation 10
            obj.shibueHConstants = zeros(3,1);
            obj.shibueHConstants(1) =  1.28746e-1;
            obj.shibueHConstants(2) = -7.31097e-1;
            obj.shibueHConstants(3) = -3.15058e2;
            
            % For Shibue (2003) equation 6 from Table 1
            obj.shibueWaterPhiConstants = zeros(6,1);
            obj.shibueWaterPhiConstants(1) = -1.42006707e-2;
            obj.shibueWaterPhiConstants(2) =  1.08369910e-2;
            obj.shibueWaterPhiConstants(3) = -1.59213160e-6;
            obj.shibueWaterPhiConstants(4) = -1.10804676e-5;
            obj.shibueWaterPhiConstants(5) = -3.14287155;
            obj.shibueWaterPhiConstants(6) =  1.06338095e-3;
        end
        function SetWagnerPrussConstants( obj )
            % For Wagner and Pruss (1993) equation 2
            obj.wagnerPrussDensityConstants = zeros(6,1);
            obj.wagnerPrussDensityConstants(1) =  1.99274064;
            obj.wagnerPrussDensityConstants(2) =  1.09965342;
            obj.wagnerPrussDensityConstants(3) = -0.510839303;
            obj.wagnerPrussDensityConstants(4) = -1.75493479;
            obj.wagnerPrussDensityConstants(5) = -45.5170352;
            obj.wagnerPrussDensityConstants(6) = -6.74694450e5;
        end
        function SetLangmuirParameters( obj )
            % For Sun and Duan (2007) equation A3 from Table A1
            obj.langmuirParameters = struct();
            
            obj.langmuirParameters.smallCageA = -24.027993;
            obj.langmuirParameters.largeCageA = -22.683049;
            
            obj.langmuirParameters.smallCageB = 3134.7529;
            obj.langmuirParameters.largeCageB = 3080.3857;
        end
        function SetCavityParameters( obj )
            obj.cavityParameters = zeros(2,1);
            obj.cavityParameters(1) = 2/46; % Munck et al 1988 - Table 1
            obj.cavityParameters(2) = 6/46; % Munck et al 1988 - Table 1
        end
        function SetDaigleParameters( obj )
            obj.daigleParameters = [ 258.4719097;    ...
                                     16.54979759;   ...
                                    -0.20037934;    ...
                                    -2.51786785;    ...
                                    -8.31210883e-2; ...
                                     2.90289187e-2; ...
                                     0.24786712;    ...
                                     5.07299816e-3; ...
                                    -1.17856658e-3; ...
                                    -8.27706806e-3 ];
        end
        %%% Utility methods
        function [ muPar , lambdaPar , xiPar ] = CalcPitzerParameters( obj , P , T )
            % Input 1-D array "c" with all 10 interaction parameters in 1->10 order
            % Output resulting summation of parameters
            for parameter = 1:3
                
                switch parameter
                    case 1
                        c = obj.pitzerConstantsMu_ChemPotential_CH4Liquid;
                    case 2
                        c = obj.pitzerConstantsLambda_2ndOrderInteraction_CH4Na;
                    case 3
                        c = obj.pitzerConstantsXi_3rdOrderInteraction_CH4NaCl;
                end
                
                switch obj.pitzerModel
                    case 1992
                        % Equation 7 in Duan et al. (1992) "The prediction of..."
                        par = c(1) ...
                            + c(2) * T ...
                            + c(3) / T ...
                            + c(4) * (T^2) ...
                            + c(5) / (680 - T) ...
                            + c(6) * P ...
                            + c(7) * P * log(T) ...
                            + c(8) * P / T ...
                            + c(9) * P / (680 - T) ...
                            + c(10) * (P^2) / T;
                    case 2006
                        % Equation 9 in Duan (2006)
                        par = c(1) ...
                            + c(2) * T ...
                            + c(3) / T ...
                            + c(4) * (T^2) ...
                            + c(5) / (T^2) ...
                            + c(6) * P ...
                            + c(7) * P * T ...
                            + c(8) * P / T ...
                            + c(9) * P / (T^2) ...
                            + c(10) * (P^2) * T;
                end
                
                switch parameter
                    case 1
                        muPar = par;
                    case 2
                        lambdaPar = par;
                    case 3
                        xiPar = par;
                end
                
            end
        end
        function [ Tr ] = CalcTr( obj , T )
            Tr = T / obj.methaneTc;
        end
        function [ Pr ] = CalcPr( obj , P )
            Pr = P / obj.methanePc;
        end
        function SetSalinityParameters( obj , salinityWtPercent )
            salinityWt = salinityWtPercent / 100; % weight % to weight fraction
            
            gramsNaClper1kgWater = salinityWt*1000 / (1 - salinityWt);
            obj.molalityNaCl = gramsNaClper1kgWater / obj.mwNaCl; % moles NaCl / 1 kg water
            
            gramsH2OAssuming1GramNaCl = (1 - salinityWt)/salinityWt;
            obj.x_NaCl = ( 1/obj.mwNaCl ) ...
                            / ( 1/obj.mwNaCl + gramsH2OAssuming1GramNaCl/obj.mwH2O );
        end
        function [ waterActivity ] = CalcWaterActivity( obj )
            
            A_phi = 0.3915; % at temperature of 25 deg C
            b = 1.2; % unit is kg1/2 / mol1/2
            
            m = obj.molalityNaCl;
            
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
            waterActivity = exp( -obj.mwH2O / 1000 * mu * m * phi);

        end
        function [ molalityCH4 ] = ConvertMoleFractionToMolality( obj , x_CH4 )
            molalityCH4 = x_CH4 * 1000 / obj.mwH2O / (1 - x_CH4);
        end
        function [ x_CH4 ] = ConvertMolalityToMoleFraction( obj , molalityCH4 )
            x_CH4 = molalityCH4 / ( molalityCH4 + 1000 / obj.mwH2O );
        end
        
        
        
        % L-G Duan 1992 - USED
        function [ bulkSolubilityLG ] = CalcBulkSolubilityLG( obj , pressure , temperature )
            % pressure is in bars
            % temperature is in K
            
            % CH4 mole fraction in the vapor phase
            y_CH4 = obj.CalcCH4GasPhaseMoleFraction( pressure , temperature );
            
            [ muPar , lambdaPar , xiPar ] = obj.CalcPitzerParameters( pressure , temperature );
            
%             model = 'PR1978';
            model = 'Duan1992';
            switch model
                case 'Duan1992'
                    Tr = obj.CalcTr( temperature );
                    Pr = obj.CalcPr( pressure );            
                    lnPhi = obj.CalcLnFugacityCoefficient( Pr , Tr );
                case 'PR1978'
                    pr1978Calculator = BCPR1978Util();
                    [ Z , A , B ] = pr1978Calculator.CalcPR78( pressure * 1e5 , temperature );
                    fugacityCoefficientPhi = pr1978Calculator.CalcFugacityCoefficient( Z , A , B );
                    lnPhi = log( fugacityCoefficientPhi );
            end

            fugacity = pressure * exp( lnPhi );
            
            bulkSolubilityLG = y_CH4 * fugacity ...
                                    * exp( -muPar - 2 * lambdaPar * obj.molalityNaCl - xiPar * obj.molalityNaCl^2 );

        end
        
        function [ y_CH4 ] = CalcCH4GasPhaseMoleFraction( obj , pressure , temperature )
            % pressure in bars
            % temperature in K
            % x_NaCl = mole fraction of NaCl in aqueous solution
            
            
            P = pressure;
            T = temperature;
            
            % bars cm^3 / mol K
            R = obj.gasConstantR;
            
            % Mole fraction of H2O in the liquid phase
            x_H2O = 1 - 2*obj.x_NaCl;
            
            % Vapor pressure of water with aqueous NaCl using Shibue (2003)
            % Units in MPa
            Psat_H2O = obj.CalcAqueousNaClVaporPressure( T );
            % Conversion from MPa to bars because rest of code is in bars
            Psat_H2O = Psat_H2O * 10;
            
            % NaCl solution Tc in K
            aqueousNaClTc = obj.CalcAqueousNaClTc();

            % NaCl solution molar volume in cm^3/mole
            if(isempty(aqueousNaClTc))
                aqueousNaClTc
            end
            liquidMolarVolume_H2O = obj.CalcLiquidMolarVolumeWater( aqueousNaClTc , T );
            
            % Fugacity coefficient of water in the gas phase
            phi_H2O = obj.CalcPhiWaterInGasPhase( pressure , temperature);
            
            y_H2O = ( x_H2O * Psat_H2O ) ...
                    / ( phi_H2O * P ) ...
                    * exp( ...
                            liquidMolarVolume_H2O * ( P - Psat_H2O ) ...
                            / ( R * T ) ...
                         );
                     
%             model = 'DaigleCode';
            model = 'DuanModel';
            switch model
                case 'DuanModel'
                    y_CH4 = 1 - y_H2O;
                case 'DaigleCode'
                    Pwsat = obj.CalcPwsat(temperature);
                    Pwsat = Pwsat * 10;
                    y_CH4 = (P - Pwsat) / P;
            end
        end
        function [ Pwsat ] = CalcPwsat( ~ , temperature )
            T = temperature;
            
            % Water properties
            Tcw = 647.14; %K
            Pcw = 22.064; %MPa
            rho_cw = 322; %kg/m^3
            % Saul & Wagner (1987) vapor pressure of water
            a1 = -7.85823;
            a2 = 1.83991;
            a3 = -11.7811;
            a4 = 22.6705;
            a5 = -15.9393;
            a6 = 1.77516;
            tau = 1-(T./Tcw);
            var1 = (a1.*tau) + (a2.*(tau.^1.5)) + (a3.*(tau.^3)) + (a4.*(tau.^3.5)) + (a5.*(tau.^4)) + (a6.*(tau.^7.5));
            % Vapor pressure in MPa
            Pwsat = Pcw.*(exp((var1.*Tcw)./T));
            deltat = 1e-5;
            taudt = 1 - ((T+deltat)./Tcw);
            var1 = (a1.*taudt) + (a2.*(taudt.^1.5)) + (a3.*(taudt.^3)) + (a4.*(taudt.^3.5)) + (a5.*(taudt.^4)) + (a6.*(taudt.^7.5));
            var2 = (var1.*Tcw)./(T+deltat);
            % Vapor pressure at T+dt in MPa
            Pwsatdt = Pcw.*exp(var2);
            dpdt = (Pwsatdt-Pwsat)./deltat;
            dpdt = dpdt.*1e6;
        end
        function [ vaporPressureAqueousNaCl ] = CalcAqueousNaClVaporPressure( obj , T )
            % T = temperature in K
            % x_NaCl = mole fraction of NaCl in aqueous solution
            
            % Shibue (2003) equation 4
            aqueousNaClPc = obj.CalcAqueousNaClPc();
            
            % Shibue (2003) equation 8
            g = obj.CalcShibueG( T );
            
            % Shibue (2003) equation 10
            h = obj.CalcShibueH();
            
            % Shibue (2003) equation 9
            lnPCalc = log( aqueousNaClPc ) + g + h;
            % Units in MPa because obj.waterPc is in MPa
            vaporPressureAqueousNaCl = exp( lnPCalc );
        end
        function [ aqueousNaClTc ] = CalcAqueousNaClTc( obj )
            % x_NaCl = mole fraction of NaCl in aqueous solution
            q = obj.shibueCriticalConstants;
            
            % Shibue (2003) equation 3
            aqueousNaClTc = obj.waterTc ...
                            + ( q(1) * obj.x_NaCl^(1/2) ) ...
                            + ( q(2) * obj.x_NaCl ) ...
                            + ( q(3) * obj.x_NaCl^2 ) ...
                            + ( q(4) * obj.x_NaCl^4 );
        end
        function [ aqueousNaClPc ] = CalcAqueousNaClPc( obj )
            % x_NaCl = mole fraction of NaCl in aqueous solution
            q = obj.shibueCriticalConstants;
            
            % Shibue (2003) equation 4
            % Units in MPa because obj.waterPc is in MPa
            aqueousNaClPc = obj.waterPc ...
                            + ( q(5)  * obj.x_NaCl   ) ...
                            + ( q(6)  * obj.x_NaCl^2 ) ...
                            + ( q(7)  * obj.x_NaCl^3 ) ...
                            + ( q(8)  * obj.x_NaCl^4 ) ...
                            + ( q(9)  * obj.x_NaCl^5 ) ...
                            + ( q(10) * obj.x_NaCl^6 );
        end
        function [ g ] = CalcShibueG( obj , T )
            % T = temperature in K
            %
            % ln(Pr) - g should range from -0.276 @ 273.15 K
            %                           to -0.290 @ 623.04 K
            %                           for x_Nacl = 0.00893 (2.84 wt%)
            
            TcOverT = obj.waterTc / T;
            TOverTc = T / obj.waterTc;
            
            g = TcOverT ...
                    * (   obj.shibueGConstants(1) * (1 - TOverTc) ...
                        + obj.shibueGConstants(2) * (1 - TOverTc)^1.5 ) ...
                + TcOverT ...
                    * (   obj.shibueGConstants(3) * (1 - TOverTc)^3 ...
                        + obj.shibueGConstants(4) * (1 - TOverTc)^3.5 ) ...
                + TcOverT ...
                    * (   obj.shibueGConstants(5) * (1 - TOverTc)^4 ...
                        + obj.shibueGConstants(6) * (1 - TOverTc)^7.5 );
        end
        function [ h ] = CalcShibueH( obj )
            % x_NaCl = mole fraction of NaCl in aqueous solution
            a = obj.shibueHConstants;
            
            h = a(2) * obj.x_NaCl ...
                / ( obj.x_NaCl + a(1)^2 ) ...
                + ( a(3) * obj.x_NaCl^2 );
        end
        function [ liquidMolarVolume_H2O ] = CalcLiquidMolarVolumeWater( obj , Tc , T )
            % Tc = critical temperature in K
            % T = temperature in K
            % x_NaCl = mole fraction of NaCl in liquid phase
            
            % NaCl solution MW in grams/mole 
            aqueousSolutionMW = obj.x_NaCl * obj.mwNaCl ...
                                + (1 - obj.x_NaCl) * obj.mwH2O;
                            
            b = obj.wagnerPrussDensityConstants;
                    
            theta = T / Tc;
            tau = 1 - theta;
            
            % Density in kg/m^3
            density = obj.wagnerPrussReferenceDensity ...
                        * ( 1 ...
                            + b(1) * tau^(1/3) ...
                            + b(2) * tau^(2/3) ...
                            + b(3) * tau^(5/3) ...
                            + b(4) * tau^(16/3) ...
                            + b(5) * tau^(43/3) ...
                            + b(6) * tau^(110/3) ...
                          );
            % conversion to g/cm^3
            density = density / 1000;
            
            % turn density into molar volume (cm^3/mol)
            liquidMolarVolume_H2O = aqueousSolutionMW / density;
            
%             aqueousSolutionMW = aqueousSolutionMW
%             density = density
        end
        function [ phi_H2O ] = CalcPhiWaterInGasPhase( obj , pressure , temperature )
            % pressure in bars
            % temperature in K
            
            P = pressure;
            T = temperature;
            
            a = obj.shibueWaterPhiConstants;
            
            phi_H2O = exp( a(1) ...
                            + a(2) * P ...
                            + a(3) * P^2 ...
                            + a(4) * P * T ...
                            + a(5) * P / T ...
                            + a(6) * (P^2) / T ...
                         );
        end
        function [ Vr , Z ] = CalcVrAndZ( obj , Pr , Tr )
            % Newton's Method iteration to calculate real Vr using
            % definition Vr = (V*Pc)/(R*Tc) // V is molar volume
            %               = (Z * T/Tc * Pc/P)
            %               = (Z * Tr / Pr)
            %
            % Only take scalar Pr and Tr, no vectorization
            
            % Perturbation for slope calculation
            eps = 1e-4; 
            
            % Initial guess of Vr using ideal gas assumption for Z = 1
            % x_n
            Vr = Tr / Pr;
            
            % Do while loop for Newton's method
            % Condition is until the difference between the two calculated
            % Z's become 0
            iteration = 0;
            doWhileFlag = true;
            while( doWhileFlag || abs(deltaZ) > 1e-10 )
                doWhileFlag = false;
                iteration = iteration + 1;
                if(iteration > 20)
                    disp('Iteration exceeded limit of 20 times for LK-Duan EOS calculation for gas bulk solubility calculation')
                    break
                end
                
                [ ZUsingReducedPVT , ZUsingEOS ] = obj.IterateDuanLKEOSForZ( Vr , Pr , Tr );
                % f(x)
                deltaZ = ZUsingReducedPVT - ZUsingEOS;
                
                [ ZUsingReducedPVTPerturbed , ZUsingEOSPerturbed ] = obj.IterateDuanLKEOSForZ( Vr + eps , Pr , Tr );
                % f(x + eps)
                deltaZPerturbed = ZUsingReducedPVTPerturbed - ZUsingEOSPerturbed;
                
                % f'(x)
                slope = (deltaZPerturbed - deltaZ) ./ eps;
                
                % x_n+1
                Vr = Vr - (deltaZ ./ slope);
            end
%             iteration
            Z = Pr * Vr / Tr;
        end
        function [ ZUsingReducedPVT , ZUsingEOS ] = IterateDuanLKEOSForZ( obj , Vr , Pr , Tr )
            % Calculates Z (compressibility factor) by using 2 definitions of
            % Z in Duan et al. (1992)'s appendix.
            %
            % Vr = reduced volume
            % Pr = P/Pc = reduced pressure
            % Tr = T/Tc = reduced temperature
            
            [ B , C , D , E , F , ~ ] = obj.CalcDuanLKEOSParameters( Tr , Vr );
            a = obj.duanLKConstants;
            
            ZUsingReducedPVT = Pr * Vr / Tr;
            
            ZUsingEOS = 1 ...
                        + B / Vr ...
                        + C / (Vr^2) ...
                        + D / (Vr^4) ...
                        + E / (Vr^5) ...
                        + F / (Vr^2) ...
                            * ( a(14) + a(15) / (Vr^2) ) ...
                            * exp( -a(15) / (Vr^2) );
        end
        function [ B , C , D , E , F , G ] = CalcDuanLKEOSParameters( obj , Tr , Vr )
            a = obj.duanLKConstants;
            
            B = a(1) ...
                + a(2) / (Tr^2) ...
                + a(3) / (Tr^3);
            
            C = a(4) ...
                + a(5) / (Tr^2) ...
                + a(6) / (Tr^3);
            
            D = a(7) ...
                + a(8) / (Tr^2) ...
                + a(9) / (Tr^3);
            
            E = a(10) ...
                + a(11) / (Tr^2) ...
                + a(12) / (Tr^3);
            
            F = a(13) / (Tr^3);
            
            G = F / (2*a(15)) ...
                * ( a(14) + 1 - ...
                        ( a(14) + 1 + a(15) / (Vr^2) ) ...
                        * exp( -a(15) / (Vr^2) ) ...
                  );
        end
        function [ lnPhi ] = CalcLnFugacityCoefficient( obj , Pr , Tr )
            
            [ Vr , Z ] = obj.CalcVrAndZ( Pr , Tr );
            [ B , C , D , E , ~ , G ] = obj.CalcDuanLKEOSParameters( Tr , Vr );
            
            % Natural log of fugacity ooefficient
            lnPhi = Z - 1 - log(Z) ...
                    + B / Vr ...
                    + C / (2 * Vr^2) ...
                    + D / (4 * Vr^4) ...
                    + E / (5 * Vr^5) ...
                    + G;
        end
        
        % L-H Henry 1999 - UNUSED/DOESN'T MATCH
        function [ bulkSolubilityLH ] = CalcBulkSolubilityLH_OLD( obj , pressure , temperature  )
            % pressure is in bars
            % temperature is in K
            
            % CH4 mole fraction in the vapor phase
            y_CH4 = obj.CalcCH4GasPhaseMoleFraction( pressure , temperature );
            
            [ muPar , lambdaPar , xiPar ] = obj.CalcPitzerParameters( pressure , temperature );

            waterActivity = obj.CalcWaterActivity();
            
            fugacity = obj.CalcFugacityLH( pressure , temperature , waterActivity );
            
            bulkSolubilityLH = y_CH4 * fugacity ...
                                    * exp( -muPar - 2 * lambdaPar * obj.molalityNaCl - xiPar * obj.molalityNaCl^2 );
        end
        function [ fugacity ] = CalcFugacityLH( obj , pressure , temperature , waterActivity )
            % pressure in bars
            % temperature in K
            
            % x_n, initial guess
            fugacity = 0.5 .* pressure;
            
            eps = 1e-4; % perturbation for slope calculation

            iteration = 0;
            doWhileFlag = true;
            while( doWhileFlag || abs( deltaChemPotential ) > 1e-10 )
                doWhileFlag = false;
                iteration = iteration + 1;
                
                if(iteration > 20)
                    disp('Iteration exceeded limit of 20 times for fugacity calculation for hydrate bulk solubility calculation')
                    fugacity = 'ERROR';
                    break
                end
                
                % f(x)
                deltaChemPotential = obj.CalcDeltaChemPotential( fugacity , pressure , temperature , waterActivity );
                
                % f(x + eps)
                deltaChemPotentialPerturbed = obj.CalcDeltaChemPotential( fugacity + eps , pressure , temperature , waterActivity );
                
                % f'(x)
                slope = (deltaChemPotentialPerturbed - deltaChemPotential) / eps;
                
                % x_n+1
                fugacity = fugacity - (deltaChemPotential / slope);
%                 deltaChemPotential
            end
        end
        function [ deltaChemPotential ] = CalcDeltaChemPotential( ~ , x , pressure , temperature , waterActivity )
            % x = fugacity in bars
            % pressure in bars
            % temperature in K
            
            P = pressure;
            T = temperature;
            R = 8.3144598; % in Pa m^3 / mol K (Gas constant)
            T_bar = (T + 273.15)/2; % Munck et al 1988 - Equation 8

            % vi is the number of cavitis of type i per molecular of water
            % 46 water, 2 small cavity and 6 large cavity per unit cell.
            v1 = 2/46; % Munck et al 1988 - Table 1
            v2 = 6/46; % Munck et al 1988 - Table 1
            
            % Munck et al 1988 - Table 6
            % Small cavity 
            A1 = 0.7228e-3 / 1.0132501; % K / atm -> converted to K / bars
            B1 = 3187; % K 
            % large cavity
            A2 = 23.35e-3 / 1.0132501; % K / atm -> converted to K / bars
            B2 = 2653; % K
            
            % Munck et al 1988 - Table 3
            delta_mu0 = 1264; % J / mol or Pa m^3 / mol
            delta_H0 = -4858; % J / mol or Pa m^3 / mol
            delta_V0 = 4.6; % cm^3 / mol
            delta_Cp = 39.16; % J / mol K or Pa m^3 / mol K 
            
%             % Sun and Duan 2007 - Table 1
%             delta_mu0 = 1202; % J / mol or Pa m^3 / mol
%             delta_H0 = 1300 - 6009.5; % J / mol or Pa m^3 / mol
            
            
            
            V_beta = 22.6;
            V_l = 18;
            
            r = 20E-9;
            gamma_hw = 0.027;
            
            % Munck et al 1988 - Equation 7, 8 &&&& Henry 1999 - Equation A4 (water activity term)
            chemPotDelta_EmptyLatticeMinusLiquidOverRT =    delta_mu0 / (R * 273.15) ...
                                                            - (delta_H0 - delta_Cp * 273.15) / R * (-1/T + 1/273.15) ...
                                                            - delta_Cp / R * log(T / 273.15) ...
                                                            + delta_V0 / (R * T_bar) * P * 0.1 ... % (* 0.1) conversion accounts for cm^3 -> m^3 and bars -> Pa
                                                            - log(waterActivity); 
                                                            
            
            
            % left=delta_mu0/(R*273.15) - ...
            %      delta_H0/R*(-1/T+1/273.15) + ...
            %      delta_V0/(R*T)*P;
            
%             V_beta * 2 * gamma_hw / r / 1e5 / (R * T)
            
            % Munck et al 1988 - Equation 4
            C1 = A1 / T * exp(B1 / T);
            C2 = A2 / T * exp(B2 / T);
            
            chemPotDelta_EmptyLatticeMinusHydrateOverRT = -( v1 * log(1 - C1 * x / (1 + C1 * x)) + v2 * log(1 - C2 * x / (1 + C2 * x)) );

            
            deltaChemPotential = chemPotDelta_EmptyLatticeMinusLiquidOverRT - chemPotDelta_EmptyLatticeMinusHydrateOverRT;
        end
        
        % L-H Duan 2006 - DOESN'T WORK
        function [ bulkSolubilityLH ] = CalcBulkSolubilityLH_NOTWORKING( obj , P , T , initialGuessSolubility )
            
            % Perturbation for slope calculation
            eps = 1e-4;
            
            % x_n
            targetSolubility = initialGuessSolubility;
            
            
            
            
            
            iteration = 0;
            doWhileFlag = true;
            while( doWhileFlag || abs( deltaMu ) > 1e-6 )
                doWhileFlag = false;
                iteration = iteration + 1;
                
%                 [ deltaMu_Water_Hydrate , deltaMu_Water_Liquid ] = obj.IterateForSolubilityLH( targetSolubility , P , T );
%                 % f(x)
%                 deltaMu = deltaMu_Water_Hydrate - deltaMu_Water_Liquid;
%                 
%                 [ deltaMu_Water_HydratePerturbed , deltaMu_Water_LiquidPerturbed ] = obj.IterateForSolubilityLH( targetSolubility + eps, P , T );
%                 % f(x + eps)
%                 deltaMuPerturbed = deltaMu_Water_HydratePerturbed - deltaMu_Water_LiquidPerturbed;
%                 
%                 % f'(x)
%                 slope = (deltaMuPerturbed - deltaMu) ./ eps;
%                 
%                 % x_n+1
%                 targetSolubility = targetSolubility - (deltaMu ./ slope);
                
                [ deltaMu_Water_Hydrate , deltaMu_Water_Liquid ] = obj.IterateForSolubilityLH( targetSolubility , P , T )
                deltaMu = deltaMu_Water_Hydrate - deltaMu_Water_Liquid
                targetSolubility = targetSolubility + 0.0001
                
                
                [ molalityCH4 ] = ConvertMoleFractionToMolality( obj , targetSolubility )
            end
            iteration
            bulkSolubilityLH = targetSolubility;
            
            
        end
        
        function [ deltaMu_Water_Hydrate , deltaMu_Water_Liquid ] = IterateForSolubilityLH_NOTWORKING( obj , targetSolubility , P , T )
            
            PsatInitialGuess = 100; % must be in bars
            Psat = obj.CalcPsat( targetSolubility , PsatInitialGuess , T );
            v = obj.cavityParameters;
            
            Tr = obj.CalcTr( T );
            Pr = obj.CalcPr( Psat );            
            lnPhi = obj.CalcLnFugacityCoefficient( Pr , Tr );
            PsatFugacity = exp(lnPhi) * P;
            
            disp('new LH fugacity')
            fugacity = obj.CalcPoyntingCorrectionFugacity( P , T , Psat , PsatFugacity )
            
            theta = obj.CalcFractionalOccupancyTheta( T , fugacity );
            summationTheta = sum(theta);
            
            
            
            deltaMu_Water_Hydrate = - obj.gasConstantR * T ...
                                    * (   v(1) * log( 1 - summationTheta ) ...
                                        + v(2) * log( 1 - summationTheta ) ...
                                      );  
            
            
            
            % refactor this later
            R = obj.gasConstantR;
            waterActivity = obj.CalcWaterActivity();
            
            delta_mu0 = 12640; % Munck et al 1988 - Table 3
            delta_H0 = -48580; % Munck et al 1988 - Table 3
            delta_V0 = 4.6; % Munck et al 1988 - Table 3
            delta_Cp = 391.6; % Munck et al 1988 - Table 3
                        
            T_bar = (T + 273.15)/2; % Munck et al 1988 - Equation 8
            
            % Munck et al 1988 - Equation 7, 8 &&&& Henry 1999 - Equation A4 (water activity term)
            deltaMu_Water_LiquidOverRT = delta_mu0 ./ (R*273.15) ... 
                                        - (delta_H0 - delta_Cp .* 273.15) ./ R .* (-1/T + 1/273.15) - delta_Cp ./ R .* (log(T) - log(273.15)) ...
                                        + delta_V0 ./ (R .* T_bar) .* P ...
                                        - log(waterActivity);
            deltaMu_Water_Liquid = deltaMu_Water_LiquidOverRT * R * T;
%             deltaMu_Water_Liquid = deltaMu_Water_LiquidOverRT;
            
            
        end
        function [ Psat ] = CalcPsat( obj , targetSolubility , PsatInitialGuess , temperature )
            
            
            % Perturbation for slope calculation
            eps = 1e-4; 
            
            % x_n
            PsatGuess = PsatInitialGuess; % must be in bars
            
            iteration = 0;
            doWhileFlag = true;
            while( doWhileFlag || abs( targetSolubility - calculatedSolubility ) > 1e-6 )
                doWhileFlag = false;
                iteration = iteration + 1;
                
                
                % f(x)
                calculatedSolubility = obj.CalcBulkSolubilityLG( PsatGuess , temperature );
                calculatedSolubility = ConvertMolalityToMoleFraction( obj , calculatedSolubility );
                
                deltaSolubility = targetSolubility - calculatedSolubility;
                
                % f(x + eps)
                calculatedSolubilityPerturbed = obj.CalcBulkSolubilityLG( PsatGuess + eps , temperature );
                calculatedSolubilityPerturbed = ConvertMolalityToMoleFraction( obj , calculatedSolubilityPerturbed );
                
                deltaSolubilityPerturbed = targetSolubility - calculatedSolubilityPerturbed;  
                
                
                % f'(x)
                slope = (deltaSolubilityPerturbed - deltaSolubility) ./ eps;
                
                % x_n+1
                PsatGuess = PsatGuess - (deltaSolubility ./ slope);
                
                
            end
%             iteration
            Psat = PsatGuess;
        end
        function [ fugacity ] = CalcPoyntingCorrectionFugacity( obj , P , T , Psat , PsatFugacity )
            % P = pressure in bars
            % T = temperature in K
            % Psat = saturation pressure in bars
            % PsatFugacity = fugacity of CH4 at saturation pressure, in bars
            fugacity = PsatFugacity * exp( ...
                                            obj.partialMolarVolumeAqueousCH4 ...
                                            * ( Psat - P ) ...
                                            / ( obj.gasConstantR * T ) ...
                                         );
        end
        function [ theta ] = CalcFractionalOccupancyTheta( obj , T , fugacity )
            % T = temperature in K
            % fugacity in bars
            
            n = 2; % number of types of cages
            
            theta = zeros(n,1);
            
            cageType = strings(n,1);
            cageType(1) = 'small';
            cageType(2) = 'large';
            
            c = zeros(n,1);
            
            for i = 1:n
                c(i) = obj.CalcLangmuirC( T , cageType(i) );
            end
            
            
            %%% testing old langmuir model
            % Small cavity
            A1 = 0.7228e-3; % K/atm % Munck et al 1988 - Table 6
            B1 = 3187; % K % Munck et al 1988 - Table 6
            
            % large cavity
            A2 = 23.35e-3; % K/atm % Munck et al 1988 - Table 6
            B2 = 2653; % K % Munck et al 1988 - Table 6
            
            c(1) = A1 ./ T .* exp(B1 ./ T); % Munck et al 1988 - Equation 4
            c(2) = A2 ./ T .* exp(B2 ./ T); % Munck et al 1988 - Equation 4 
            
            
            summationCFugacity = sum( c .* fugacity );
            
            for i = 1:n
                theta(i) = c(i) * fugacity ...
                            / ( 1 + summationCFugacity );            
            end
            
        end
        function [ C ] = CalcLangmuirC( obj , T , cageType )
            % T = temperature in K
            
            % Sun and Duan (2007) equation A3
            switch cageType
                case 'small'
                    C = exp(   obj.langmuirParameters.smallCageA ...
                             + obj.langmuirParameters.smallCageB / T ...
                           );
                case 'large'
                    C = exp(   obj.langmuirParameters.largeCageA ...
                             + obj.langmuirParameters.largeCageB / T ...
                           );
            end
        end
        
        % L-H Daigle (Davie and Buffett) - USED
        function [ bulkSolubilityLH ] = CalcBulkSolubilityLH( obj , gasBulkSolubilityT3P , temperature , T3PArray )
            bulkSolubilityLH = gasBulkSolubilityT3P .* exp( (temperature - T3PArray) ./ obj.davieCorrectionFactor );
        end
        function [ T3P ] = CalcT3P( obj , pressure , salinityWtPercent )
            % can make a = obj.daigleParameters and set daigleParameters
            % pressure in MPa
            % depth in meters
            % salinity in weight %
            % pressureGradient in MPa/m
            
            a = obj.daigleParameters;
            P = pressure;
            salt = salinityWtPercent;
            
            T3P =   a(1)                            ...
                    + ( a(2) .* log(P) )            ...
                    + ( a(3) .* salt )              ...
                    + ( a(4) .* log(P).^2 )         ...
                    + ( a(5) .* salt.^2 )           ...
                    + ( a(6) .* log(P) .* salt )    ...
                    + ( a(7) .* log(P).^3 )         ...
                    + ( a(8) .* salt.^3 )           ...
                    + ( a(9) .* log(P) .* salt.^2 ) ...
                    + ( a(10) .* salt .* log(P).^2 );
            
        end
        % UNUSED
        function [ dfdz ] = CalcDfdz( obj , P0 , salinityWtPercent , depth , tempGradient )
            a = obj.daigleParameters;
            P = P0 + 1e-2 * depth;
            dfdz = tempGradient - 1e-2/P ...
                                * ( a(2) ...
                                    + 2*a(4)*log(P) ...
                                    + a(6)*salinityWtPercent ...
                                    + 3*a(7)*(log(P))^2 ...
                                    + a(9)*(salinityWtPercent)^2 ...
                                    + 2*a(10)*salinityWtPercent*log(P) ...
                                  );
        end
        function [ T3PArray ] = CalcT3PArray( ~ , pressure )
            P = pressure;
            T3PArray = 1./(-3.6505e-6.*(log(P)).^3+3.704e-5.*(log(P)).^2 - 0.00021948.*log(P) + 0.0038494);
        end
        function [ depth3PBulkEQL , temp3PBulkEQL ] = CalcDepth3PBulkEQL( obj , P0 , T0 , salinityWtPercent , temperatureGradient )
            P0 = P0 / 1e6; % change Pa to MPa
            G = temperatureGradient / 1000; % change C/km to C/m
            d = 100; % initial guess of 3P bulk EQL depth in meters
            
            T3P = obj.CalcT3P( P0 , salinityWtPercent , d );
            fz = T0 + (G * d) - T3P;
            norm2 = abs(fz);
            norm1 = 1e6;
            
            iteration = 0;
            while( norm1 > 1e-8 || norm2 > 1e-8 )
                iteration = iteration + 1;
                dfdz = obj.CalcDfdz( P0 , salinityWtPercent , d , G );
                dNew = d - (fz/dfdz);
                norm1 = abs(d - dNew);
                T3P = obj.CalcT3P( P0 , salinityWtPercent , dNew );
                fz = T0 + (G * dNew) - T3P;
                norm2 = abs(fz);
                d = dNew;
            end
            
            depth3PBulkEQL = d;
            temp3PBulkEQL = T0 + (G * depth3PBulkEQL);
        end

        %%% Main methods
        function [ bulkSolubilityLG , bulkSolubilityLH , T3PArray ] = CalcBulkSolubilities( obj , pressurePa , temperature , salinityWtPercent )
            % pressure in Pa
            % temperature in K
            
            n = numel(pressurePa);
            bulkSolubilityLG = zeros(n,1);
            gasBulkSolubilityT3P = zeros(n,1);
            
            pressureBars = pressurePa .* 1e-5; % P converted from Pa to bars
            obj.SetSalinityParameters( salinityWtPercent );
            pressureMPa = pressureBars ./ 10;
            
            T3PArray = obj.CalcT3P( pressureMPa , salinityWtPercent );
            
            for i = 1:n
                bulkSolubilityLG(i) = obj.CalcBulkSolubilityLG( pressureBars(i) , temperature(i) );
                gasBulkSolubilityT3P(i) = obj.CalcBulkSolubilityLG( pressureBars(i) , T3PArray(i) );
            end
            bulkSolubilityLH = obj.CalcBulkSolubilityLH( gasBulkSolubilityT3P, temperature, T3PArray );
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
    
        
        
        
        
        
end




























































































