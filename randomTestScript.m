clc
clear variables
close all


sg = linspace(0, 0.5, 100);
sh = linspace(0, 0.5, 100);

[sgMesh, shMesh] = meshgrid(sg, sh);

competitionFractionOfGas = 0.75;
competitionFractionOfHydrate = 1 - competitionFractionOfGas;

nSg = numel(sg);
nSh = numel(sh);

adjustedSg = zeros(nSg, nSh);
adjustedSh = zeros(nSg, nSh);


for iSg = 1:nSg
    for iSh = 1:nSh
        % Above equal pore size invasion point
        if( sh(iSh)/(sh(iSh) + sg(iSg)) > competitionFractionOfHydrate)
            adjustedSg(iSg, iSh) = sg(iSg) / competitionFractionOfGas;
            adjustedSh(iSg, iSh) = sh(iSh) + sg(iSg);
        else
            adjustedSg(iSg, iSh) = sg(iSg) + sh(iSh);
            adjustedSh(iSg, iSh) = sh(iSh) / competitionFractionOfHydrate;
        end
    end
end

figure
surf(sgMesh, shMesh, adjustedSg)
figure
surf(sgMesh, shMesh, adjustedSh)





