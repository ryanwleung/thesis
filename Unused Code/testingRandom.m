clc
clear variables
close all

% asdf = CSolubilityUtil();
% asdf.SetSalinityParameters( 3.5 );
% solubility = 0.1;
% x_CH4 = asdf.ConvertMolalityToMoleFraction( solubility )
% asdf.CalcPsat(x_CH4,200,277)


asdf = CSolubilityUtil();
% asdf.SetSalinityParameters( 14.91679 );
asdf.SetSalinityParameters( 3.5 );
% asdf.x_NaCl
% asdf.molalityNaCl


% pressure = (0:2:800)';
pressure = 70;
temperature = 298;
% temperature = 511;


n = numel(pressure);
y_H2O_3m = zeros(n,1);
for i = 1:n
    y_H2O_3m(i) = 1 - asdf.CalcCH4GasPhaseMoleFraction( pressure(i) , temperature );
end
figure
plot(pressure,y_H2O_3m,'--')
    
asdf.SetSalinityParameters( 0 );
y_H2O_0m = zeros(n,1);
for i = 1:n
    y_H2O_0m(i) = 1 - asdf.CalcCH4GasPhaseMoleFraction( pressure(i) , temperature );
end
hold on
plot(pressure,y_H2O_0m)

% axis([0 400 -0.001 0.006])
axis([0 800 0 0.5])
    
    
    
    
    