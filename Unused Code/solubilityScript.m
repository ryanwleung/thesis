clc
close all
clear variables


depth = obj.DataTable.depth;


pressure = obj.CalcPressure( depth );
temperature = obj.CalcTemperature( depth );
gasDensity = BCFormation.CalcGasDensityArray( pressure , temperature );

obj.DataTable.gasBulkSolubilityInterpolated = interp1( obj.Bulk.Depth , obj.Bulk.Solubility , depth );



solubilityCalculator = BCSolubilityUtil();

[ gasBulkSolubility , hydrateBulkSolubility , T3PArray ] = ...
    solubilityCalculator.CalcBulkSolubilities( pressure , temperature , obj.salinityWtPercent );
        
% % BR adjustment to match Liu bulk sol results
% gasBulkSolubility = gasBulkSolubility - .002;
% hydrateBulkSolubility = hydrateBulkSolubility + .0024;



interpolatedBulkSolLogical = ~isnan(obj.DataTable.gasBulkSolubilityInterpolated);
tempInterpBulkSol = obj.DataTable.gasBulkSolubilityInterpolated(interpolatedBulkSolLogical);
tempInterpDepth = obj.DataTable.depth(interpolatedBulkSolLogical);

figure(1)
switch class(obj)
    case 'DCHydrateRidge'
        tempInterpDepth = tempInterpDepth + 790;
        depth = depth + 790;
    case 'DCBlakeRidge'
        tempInterpDepth = tempInterpDepth + 2780;
        depth = depth + 2780;
end


plot( tempInterpBulkSol , tempInterpDepth , 'k' , 'LineWidth' , 3 )
hold on
plot( gasBulkSolubility , depth , 'r-' )
hold on
plot( hydrateBulkSolubility , depth , 'g-' )
hold on
xlabel('Solubility')
ylabel('Depth (mbsf)')
set(gca,'YDir','Reverse')
switch class(obj)
    case 'DCHydrateRidge'
        title('Hydrate Ridge')
        axis([0.08 0.125 850 950])
    case 'DCBlakeRidge'
        title('Blake Ridge')
        axis([0.155 0.205 3200 3300])
end


switch class(obj)
    case 'DCHydrateRidge'
        load('Hydrate Ridge Data\Solubility Plots\BulkSolHR.mat')
        plot( BulkSolHR(:,1) , BulkSolHR(:,2) )    
    case 'DCBlakeRidge'
        load('Blake Ridge Data\Solubility Plots\BulkSolBR.mat')
        plot( BulkSolBR(:,1) , BulkSolBR(:,2) )
end





% return

ch4Quantity = 40;

% for ch4Quantity = [6 12 18 24 30 36]
    



[ gasMaxSolubility , gasSaturation2P ] = obj.CalcMaxSolLG( ch4Quantity , pressure , gasDensity , gasBulkSolubility );
[ hydrateMaxSolubility , hydrateSaturation2P ] = obj.CalcMaxSolLH( ch4Quantity , temperature , hydrateBulkSolubility );
% gasMaxSolubility = gasMaxSolubility + 0.0062;

plot( gasMaxSolubility , depth , 'r--' )
hold on

plot( hydrateMaxSolubility , depth , 'g--' )
hold on

gasMinSolubility = BCFormation.CalcSolubilityLG( gasBulkSolubility , obj.CalcPcgw(0) * 1e6 , pressure );
hydrateMinSolubility = BCFormation.CalcSolubilityLH( hydrateBulkSolubility , obj.CalcPchw(0) * 1e6 , temperature );

plot( gasMinSolubility , depth , 'r-.' )
hold on

plot( hydrateMinSolubility , depth , 'g-.' )
hold on


% return


top3PIndex = BCFormation.GetTop3PIndex( gasMinSolubility , hydrateMaxSolubility );
bottom3PIndex = BCFormation.GetBottom3PIndex( gasMaxSolubility , hydrateMinSolubility );
thickness = BCFormation.GetThickness3PZone( depth , top3PIndex , bottom3PIndex );


sg = gasSaturation2P;
sg( 1 : bottom3PIndex - 1 ) = 0;

sh = hydrateSaturation2P;
sh( top3PIndex + 1 : end ) = 0;

sol = hydrateMaxSolubility;
sol( bottom3PIndex : end ) = gasMaxSolubility( bottom3PIndex : end );

indexArrayOf3PZone = top3PIndex + 1 : bottom3PIndex - 1;


[ sg3P , sh3P , sol3P ] = obj.Calc3P( ch4Quantity , indexArrayOf3PZone , ...
                                    pressure , temperature , gasDensity , ...
                                    gasBulkSolubility , hydrateBulkSolubility , hydrateMaxSolubility(top3PIndex) );

sg(indexArrayOf3PZone) = sg3P;
sh(indexArrayOf3PZone) = sh3P;
sol(indexArrayOf3PZone) = sol3P;

hold on
plot( sol , depth , 'b' )



figure
plot(sg,depth,'r')
hold on
plot(sh,depth,'g')
xlabel('Saturation')
ylabel('Depth (mbsf)')
set(gca,'YDir','Reverse')
switch class(obj)
    case 'DCHydrateRidge'
        title('Hydrate Ridge')
        axis([0 0.5 850 950])
    case 'DCBlakeRidge'
        title('Blake Ridge')
        axis([0 0.33 3200 3300])
end

% end

thickness


% 
% figure
% 
% plot(gasSaturation2P,depth)
% hold on
% plot(gasSaturationBulk2P,depth)
% set(gca,'YDir','Reverse')
% axis([0 0.33 3200 3300])




























