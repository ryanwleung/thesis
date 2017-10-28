function [ Vr , B , C , D , E , G ] = calcVr( Pr , Tr )
% Newton's Method iteration to calculate real Vr using definition Vr = (V*Pc)/(R*Tc)
% Only take scalar Pr and Tr, no vectorization

eps = 1e-4; % perturbation for slope calculation

% Initial guess of Vr using IG assumption for Z = 1
Vr = Tr/Pr; 
[ deltaZ ] = calcDuanEOS( Vr , Pr , Tr );

% iteration = 0;
while abs(deltaZ) > 1e-6
%     iteration = iteration + 1
    dDeltaZ_dVr = (calcDuanEOS( Vr + eps , Pr , Tr ) - deltaZ)/eps;
    Vr_new = Vr - (deltaZ / dDeltaZ_dVr);
    
    Vr = Vr_new;
    [ deltaZ , B , C , D , E , G ] = calcDuanEOS( Vr , Pr , Tr );
    
end