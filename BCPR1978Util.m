classdef BCPR1978Util < handle
    % Peng-Robinson 1978 EOS
    
    properties (Constant)
        % Methane attributes
        acentricFactorCH4 = .008;
        mwCH4 = 16.043; % g/mol
        tcCH4 = 190.6; % K
        pcCH4 = 4.6 * 1e6; % Pa
        gasConstantR = 8.314; % in Pa m^3 / mol K (Gas constant)
    end
    
    methods
        function [ gasDensity ] = CalcGasDensity( obj , pressure , temperature )
            % pressure in Pa
            % temperature in K
            
            % Assumes the correct z compressibility factor is z1
            [ z1 ] = obj.CalcPR78( pressure , temperature );
            z1 = 1;
            gasDensity = (pressure * obj.mwCH4)/(z1 * obj.gasConstantR * temperature) / 1000; ... in kg/m^3
        end
        function [ Z , A , B ] = CalcPR78( obj , P , T )
            w = obj.acentricFactorCH4;
            Tc = obj.tcCH4;
            Pc = obj.pcCH4;
            R = obj.gasConstantR;
            
            
            if( w <= .49)
                K = .37464 + 1.54226*w - .26992*w^2;
            else
                K = .37964 + w*(1.48503 + w*(-.164423 + .01666*w));
            end
            alpha = (1 + K*(1-(T/Tc)^(1/2)))^2;
            a = .457236 * R^2 * Tc^2 * alpha / Pc;
            b = .0778 * R * Tc / Pc;
            da_dt = -.45724 * R^2 * Tc^2 / Pc * K * (alpha / (T * Tc))^(1/2);
            
            A = a*P/(R*T)^2;
            B = b*P/R/T;
            a0 = -A*B + B^2 + B^3;
            a1 = A - 3*B^2 - 2*B;
            a2 = -1 + B;
            
            Q = (3*a1 - a2^2)/9;
            R = (9*a2*a1 - 27*a0 - 2*a2^3)/54;
            D = Q^3 + R^2;
            S = (R + D^(1/2))^(1/3);
            T = (R - D^(1/2))^(1/3);
            
            if(D<0)
                theta = acos(R/((-(Q^3))^(1/2)));
                
                z1 = 2*((-Q)^(1/2))*cos(theta/3) - (1/3)*a2;
                z2 = 2*((-Q)^(1/2))*cos((theta + 2*pi)/3) - (1/3)*a2;
                z3 = 2*((-Q)^(1/2))*cos((theta + 4*pi)/3) - (1/3)*a2;
            else
                z1 = (-1/3)*a2 + (S + T);
                z2 = (-1/3)*a2 - (1/2)*(S + T) + (1/2)*1i*3^(1/2)*(S - T);
                z3 = (-1/3)*a2 - (1/2)*(S + T) - (1/2)*1i*3^(1/2)*(S - T);
            end
            
            Z = z1;
        end
        function [ fugacityCoefficientPhi ] = CalcFugacityCoefficient( ~ , Z , A , B )
            % Compute fugacity
            fugacityCoefficientPhi = exp( Z - 1 - log(Z - B) ...
                                            - ( A / (sqrt(8) * B) ) ...
                                                * log(  (Z + (1 + sqrt(2)) * B) ...
                                                       /(Z + (1 - sqrt(2)) * B) ...
                                                     ) ...
                                         );
        end
    end
    
end


















































