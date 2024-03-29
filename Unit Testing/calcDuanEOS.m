function [ deltaZ , B , C , D , E , G ] = calcDuanEOS( Vr , Pr , Tr )
% Calculates the delta Z (compressibility factor) by the 2 definitions of
% Z in Duan et al. (1992)'s appendix. The correct Z and Vr is found when
% delta Z is returned as ~0

a1 =    8.72553928E-2;
a2 =    -7.52599476E-1;
a3 =    3.75419887E-1;
a4 =    1.07291342E-2;
a5 =    5.49626360E-3;
a6 =	-1.84772802E-2;
a7 =	3.18993183E-4;
a8 =	2.11079375E-4;
a9 =	2.01682801E-5;
a10 =	-1.65606189E-5;
a11 =	1.19614546E-4;
a12 =	-1.08087289E-4;
alpha =	4.48262295E-2;
beta =	7.53970000E-1;
gamma =	7.7167000E-2;

B = a1 + a2/Tr^2 + a3/Tr^3;
C = a4 + a5/Tr^2 + a6/Tr^3;
D = a7 + a8/Tr^2 + a9/Tr^3;
E = a10 + a11/Tr^2 + a12/Tr^3;
F = alpha/Tr^3;
G = F/(2*gamma)*(beta + 1 - (beta + 1 + gamma/Vr^2)*exp(-gamma/Vr^2));

deltaZ = (Pr*Vr/Tr) - ( 1 + B/Vr + C/Vr^2 + D/Vr^4 + E/Vr^5 + F/Vr^2*(beta + gamma/Vr^2)*exp(-gamma/Vr^2) );
