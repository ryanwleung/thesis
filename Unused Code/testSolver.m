n = numel(indexArrayOf3PZone);
sg3P = zeros(n,1);
sh3P = zeros(n,1);
adjustedSol = zeros(n,1);

% Start of Newton's method

% Perturbation for slope calculation
eps = 1e-4;
%             eps = -1e-4;

% Initial guess of solubility
% solubility = hydrateMaxSolubilityAtTop;
solubility = .105;

% Separate indexing for saving calculated sg, sh, solubility
%             for i = 1:n
for i = n:-1:1

    i3P = indexArrayOf3PZone(i);

    % Get previous Sg for inital guess of Sg

%                 if( i == 1)
%                     sg = 0;
%                 else
%                     sg = sg3P(i - 1);
%                 end

    if( i == n )
        % hydrate ridge bottom sg
        sg = 0.495615110897242;

        % blake ridge bottom sg
%                     sg = 0.157738524846419;
    else
        sg = sg3P(i + 1);
    end

    fzero('calc3PNewtonIteration', sg , optimset('Display','on') , obj , ch4Quantity , ...
                                                                                solubility , ...
                                                                                pressure , temperature , gasDensity , ...
                                                                                gasBulkSolubility , hydrateBulkSolubility )

    

    

    iteration
%                 sg3P(i3P) = sg;
%                 sh3P(i3P) = sh;
%                 adjustedSol(i3P) = solubility;
%                 i3P = i3P + 1;
    sg3P(i) = sg;
    sh3P(i) = sh;
    adjustedSol(i) = solubility;
end
            
            
           