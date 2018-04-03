clc
clear variables
close all


radius = linspace(0, 100, 1001);




peak = 20;
lower = 5;
upper = 80;
sigma1 = .6;
mu1 = log(peak) + sigma1^2;
% mean = exp(mu + sigma^2 / 2);
logncdf(upper, mu1, sigma1) - logncdf(lower, mu1, sigma1)
% variance1 = (exp(sigma1^2) - 1) * exp(2 * mu1 + sigma1^2);


pdf1 = lognpdf(radius, mu1, sigma1);
cdf1 = logncdf(radius, mu1, sigma1);

peak = 7;
lower = 0.01;
upper = 80;
sigma2 = .93;
mu2 = log(peak) + sigma2^2;
% mean = exp(mu + sigma^2 / 2);
logncdf(upper, mu2, sigma2) - logncdf(lower, mu2, sigma2)
pdf2 = lognpdf(radius, mu2, sigma2);
cdf2 = logncdf(radius, mu2, sigma2);


% logninv(0.05, mu, sigma)


figure
hold on
plot(radius, pdf1)
plot(radius, pdf2)
xlabel('Pore radius (microns)')
ylabel('PDF')
legend(strcat('mu=', string(mu1), '   ', 'sigma=', string(sigma1)), strcat('mu=', string(mu2), 'sigma=', string(sigma2)))


figure
hold on
plot(radius, cdf1)
plot(radius, cdf2)
xlabel('Pore radius (microns)')
ylabel('CDF')
legend(strcat('mu=', string(mu1), '   ', 'sigma=', string(sigma1)), strcat('mu=', string(mu2), 'sigma=', string(sigma2)))

figure
hold on
plot(1 - cdf1, radius)
plot(1 - cdf2, radius)
xlabel('Sw')
ylabel('Pore radius (microns)')
legend(strcat('mu=', string(mu1), '   ', 'sigma=', string(sigma1)), strcat('mu=', string(mu2), 'sigma=', string(sigma2)))




% obj = DCTheoreticalFormation(1000, 0.4);
% 
% n = 101;
% sg = linspace(0, 1, n)';
% sgIterated = zeros(n, 1);
% solLGIterated = zeros(n, 1);
% 
% for i = 1:n
%     [solLGIterated(i), sgIterated(i)] = obj.CalcMaxSolLGIteration(sg(i), 0.131725784500610, 40, 1.003982286827520e+07, 70.102118792622761);
% end
% 
% figure
% hold on
% plot(sg, sgIterated)
% plot([0, 1], [0, 1])
% plot(sg, sgIterated - sg)
% axis([0, 1, 0, 1])

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





