function [ TableExport , StorageExport , obj ] = mainFunction( formation )
    
    switch formation
        case 'Hydrate Ridge'
            obj = HydrateRidge();
        case 'Blake Ridge'
            obj = BlakeRidge();
    end




    nDepth = height(obj.DataTable);
    zerosArray = zeros(nDepth,1);
    
    [ obj.DataTable.pressure ] = obj.calcPressure();
    [ obj.DataTable.temperature ] = obj.calcTemperature();
    [ obj.DataTable.gasDensity ] = Formation.calcGasDensityArray( obj.DataTable.pressure , obj.DataTable.temperature );
    
    
    obj.DataTable.gasSolubility = interp1( obj.Bulk.Depth , obj.Bulk.Solubility , obj.DataTable.depth );

    % % WORK ON THIS CODE NEXT
    % Depth.C_L_H_Top = 100;    % m
    % Depth.C_L_H_Bottom = 540; % m, from Guerin 
    % % Depth.C_L_H_Bottom = 520; % m
    % [ C_eq , T_array ] = Calc_CH4_C_L_H( Depth.C_L_H_Bottom , Depth.C_L_H_Top , Depth.Seafloor , Con.Water_Density , Bool );
    % 
    % Data.log(:,11) = Formation.arrayInsert( size(Data.log,1) , find(Data.log(:,1) == Depth.C_L_H_Top) - Data.depth_top , C_eq );
    % % WORK ON THIS CODE NEXT

            % Solubility calculations - skipping for now

            % NaCl_MW = 58.44277; % g/mol
            % NaCl_molar = salinityWt / 100 * 1000 / NaCl_MW; % in mol NaCl/L H2O

            % CH4_L_G = Calc_CH4_C_L_G( Depth_Pressure , Depth_Temp , NaCl_molar );
            % CH4_L_G = Calc_CH4_C_L_G( Depth_Pressure , Depth_Temp , .6 );






    % Loops from 6->40 g/dm^3 (=index-1)
    methaneQuantityArray = 6:40;

    zerosMatrix = zeros(obj.Data.depth_interval,length(methaneQuantityArray));
    
    obj.Storage.sg2P = zerosMatrix;
    obj.Storage.pcgw2P = zerosMatrix;
    obj.Storage.ratio2P = zerosMatrix;

    obj.Storage.sg3P = zerosMatrix;
    obj.Storage.pcgw3P = zerosMatrix;
    obj.Storage.ratio3P = zerosMatrix;
    


    for iMethaneQuantity = 1:length(methaneQuantityArray)
        methaneQuantity = methaneQuantityArray(iMethaneQuantity);


        obj.DataTable.hydrateSaturation = Formation.arrayInsert( zerosArray , obj.Depth.saturationTop , obj.Saturation.Hydrate_Data(:,methaneQuantity + 1) );
        obj.DataTable.gasSaturation = Formation.arrayInsert( zerosArray , obj.Depth.Graph_Top , obj.Saturation.Gas_Data(:,methaneQuantity + 1) );
        obj.DataTable.waterSaturation = 1 - obj.DataTable.hydrateSaturation - obj.DataTable.gasSaturation;


        obj.DataTable.bulkDensity = obj.estimateBulkDensity();
        if(sum( ismember(properties(obj) , 'LDT' )) == 1)
            obj.DataTable.bulkDensity = Formation.arrayInsert( obj.DataTable.bulkDensity , obj.LDT.Fixed_Depth_TOP , obj.LDT.Fixed_RHOB );
        end


        [ obj.DataTable.rockStrength ] = obj.calcRockStrength();
        
        
        

        
        % Sg = calculated from mass balance // 2-Phase - Bulk Equilibrium Case 
        obj.DataTable.sg2P = Formation.calcSg2P( methaneQuantity , obj.DataTable.gasSolubility , obj.waterDensity, obj.DataTable.gasDensity );

        
        
        [ top3PIndex ] = obj.getTop3PIndex();
        [ mid3PIndex ] = obj.getMid3PIndex();
        % Sg' = 2*S_g // 3-Phase - Liu and Flemings
        obj.DataTable.sg3P = zerosArray;
        obj.DataTable.sg3P( top3PIndex : mid3PIndex - 1 ) = obj.DataTable.gasSaturation( top3PIndex : mid3PIndex - 1 ) .* 2;
        obj.DataTable.sg3P( mid3PIndex : end) = obj.DataTable.gasSaturation( mid3PIndex : end) + obj.DataTable.hydrateSaturation( mid3PIndex : end);

        [ obj.DataTable.pcgw2P ] = calcPcgw( obj , obj.DataTable.sg2P );
        [ obj.DataTable.pcgw3P ] = calcPcgw( obj , obj.DataTable.sg3P );

        obj.DataTable.ratio2P = obj.DataTable.pcgw2P ./ obj.DataTable.rockStrength;
        obj.DataTable.ratio3P = obj.DataTable.pcgw3P ./ obj.DataTable.rockStrength;


        % 14 ---> Sg' = Sg + Sh //Theoretical Max Case
        % Data.log(Depth.Top_3P_Index:Depth.Mid_3P_Index - 1,14) = Data.log(Depth.Top_3P_Index:Depth.Mid_3P_Index - 1,7) + Data.log(Depth.Top_3P_Index:Depth.Mid_3P_Index - 1,6);

        % 15 ---> Sg' = S_g + [(z - z_top)/delta_z]*S_h //Linear scaling
        % Data.log(Depth.Top_3P_Index:Depth.Mid_3P_Index - 1,15) = Data.log(Depth.Top_3P_Index:Depth.Mid_3P_Index - 1,7) + ((Data.log(Depth.Top_3P_Index:Depth.Mid_3P_Index - 1,1) - Data.log(Depth.Top_3P_Index,1))./Con.Delta_Depth) .* Data.log(Depth.Top_3P_Index:Depth.Mid_3P_Index - 1,6);


        
        obj.Storage.sg2P(:,iMethaneQuantity) = obj.DataTable.sg2P;
        obj.Storage.pcgw2P(:,iMethaneQuantity) = obj.DataTable.pcgw2P;
        obj.Storage.ratio2P(:,iMethaneQuantity) = obj.DataTable.ratio2P;

        obj.Storage.sg3P(:,iMethaneQuantity) = obj.DataTable.gasSaturation;
        obj.Storage.pcgw3P(:,iMethaneQuantity) = obj.DataTable.pcgw3P;
        obj.Storage.ratio3P(:,iMethaneQuantity) = obj.DataTable.ratio3P; 
        
        TableExport = obj.DataTable;
        StorageExport = obj.Storage;
        
    end
    
%     obj.generateResultPlots()



end