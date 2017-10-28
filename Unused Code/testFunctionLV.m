clc
clear variables


obj = BCSolubilityUtil();
obj.SetSalinityParameters( 3.5 );
P0 = 7900000 / 1e6;
T0 = 4 + 273.15;
G = 59/1000;
P = 1;

depth = 1:1:200;

mNa = obj.molalityNaCl;
% iterate to find base of MHSZ
T_try = 0;
P_try = 0;
rhog = 0;
Salt = 100*((mNa*0.05844247)/(1+(mNa*0.05844247)));
cm_triplept = 0;
a1 = 258.4719097;
a2 = 16.54979759;
a3 = -0.20037934;
a4 = -2.51786785;
a5 = -8.31210883e-2;
a6 = 2.90289187e-2;
a7 = 0.24786712;
a8 = 5.07299816e-3;
a9 = -1.17856658e-3;
a10 = -8.27706806e-3;
d = 100;
T3P = a1 + (a2*log(P0+1e-2*d)) + (a3*Salt) + (a4*(log(P0+1e-2*d))^2) + (a5*(Salt^2)) + (a6*(log(P0+1e-2*d))*Salt) + (a7*(log(P0+1e-2*d))^3) + (a8*(Salt)^3) + (a9*log(P0+1e-2*d)*(Salt)^2) + (a10*Salt*(log(P0+1e-2*d))^2);
T3P
fz = T0 + (G*d) - T3P;
norm2 = abs(fz);
norm1 = 1e6;
while norm1>1e-10 && norm2>1e-10
    d
    dfdz = G - (1e-2/(P0+1e-2*d))*(a2 + 2*a4*log(P0+1e-2*d) + a6*Salt + 3*a7*(log(P0+1e-2*d))^2 +  a9*(Salt)^2 + 2*a10*Salt*log(P0+1e-2*d));
    dnew = d - (fz/dfdz);
    norm1 = abs(d-dnew);
    T3P = a1 + a2*log(P0+1e-2*dnew) + a3*Salt + a4*(log(P0+1e-2*dnew))^2 + a5*(Salt^2) + a6*(log(P0+1e-2*dnew))*Salt+a7*(log(P0+1e-2*dnew))^3 + a8*(Salt)^3 + a9*log(P0+1e-2*dnew)*(Salt)^2 + a10*Salt*(log(P0+1e-2*dnew))^2;
    fz = T0 + (G*dnew) - T3P;
    norm2 = abs(fz);
    d = dnew;
end
% set base of MHSZ
Lt = dnew;
% compute temp at base of MHSZ
Teq = T0 + (G*dnew);

% compute L-V solubility profile
liquid_vapor_solubility;
LV_sol = sol;
clear T
T3P = 1./(-3.6505e-6.*(log(P)).^3+3.704e-5.*(log(P)).^2 - 0.00021948.*log(P) + 0.0038494);
T=T3P;
liquid_vapor_solubility;
LH_sol = sol;
clear T
T = T0 + (depth.*G);
% inside the MHSZ, compute solubility from Buffett's expression
% below MHSZ, use L-V solubility
for i=1:length(depth)
    if T(i)>Teq
        cm_eq(i)=LV_sol(i);
    else
        cm_eq(i)=LH_sol(i)*exp((T(i)-T3P(i))/15.3);
    end
end
% smooth the solubility curve at the base of MHSZ
count = length(depth) - ( round( (Lt / (max(depth)) ) * (length(depth)) ) );
corr=abs(cm_eq(count)-cm_eq(count+1));
for i=1:length(depth)
    if T(i)<Teq
        cm_eq(i)=cm_eq(i)+corr;
    else
        cm_eq(i)=cm_eq(i);
    end
end
% convert to mass fraction
for i=1:length(depth)
    if T(i)<Teq
        cm_hyd(i) = (16.*cm_eq(i))./(18 - (2.*cm_eq(i)));
    else
        cm_hyd(i) = (16.*max(cm_eq))./(18 - (2.*max(cm_eq)));
    end
end
% compute gas density from ideal gas law
rhog = 1000*(P.*16)./(8.314.*T);
% compute gas viscosity from Lennard-Jones parameters
mug = (((16.*T).^0.5)./(((2.44.*((T3P./(P.*9.86923267)).^(1/3))).^2).*1.401)).*2.6693e-6;