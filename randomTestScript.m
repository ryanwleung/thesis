clc
clear variables
close all





obj = DCTheoreticalFormation(1000, 0.4);

n = 101;
sg = linspace(0, 1, n)';
sgIterated = zeros(n, 1);
solLGIterated = zeros(n, 1);

for i = 1:n
    [solLGIterated(i), sgIterated(i)] = obj.CalcMaxSolLGIteration(sg(i), 0.131725784500610, 40, 1.003982286827520e+07, 70.102118792622761);
end

figure
hold on
plot(sg, sgIterated)
plot([0, 1], [0, 1])
plot(sg, sgIterated - sg)
axis([0, 1, 0, 1])

% sg = linspace(0, 0.5, 100);
% sh = linspace(0, 0.5, 100);
% 
% [sgMesh, shMesh] = meshgrid(sg, sh);
% 
% competitionFractionOfGas = 0.75;
% competitionFractionOfHydrate = 1 - competitionFractionOfGas;
% 
% nSg = numel(sg);
% nSh = numel(sh);
% 
% adjustedSg = zeros(nSg, nSh);
% adjustedSh = zeros(nSg, nSh);
% 
% 
% for iSg = 1:nSg
%     for iSh = 1:nSh
%         % Above equal pore size invasion point
%         if( sh(iSh)/(sh(iSh) + sg(iSg)) > competitionFractionOfHydrate)
%             adjustedSg(iSg, iSh) = sg(iSg) / competitionFractionOfGas;
%             adjustedSh(iSg, iSh) = sh(iSh) + sg(iSg);
%         else
%             adjustedSg(iSg, iSh) = sg(iSg) + sh(iSh);
%             adjustedSh(iSg, iSh) = sh(iSh) / competitionFractionOfHydrate;
%         end
%     end
% end
% 
% figure
% surf(sgMesh, shMesh, adjustedSg)
% figure
% surf(sgMesh, shMesh, adjustedSh)





