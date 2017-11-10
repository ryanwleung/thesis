classdef DCSeismicAnalysisBR < DCBlakeRidge
    properties
        graphTop
        graphBottom
        
        SaturationLF
        Dickens
    end
    properties (Constant)
        referenceTemperature = 18; % deg C
        liuTop = 3200; % mbsl
        liuBottom = 3300; % mbsl
        saturationTop = 200; % mbsf
        saturationBottom = 650; % mbsf
        
        clayK = 21.2*10^9; % Pa
    end
    methods
        function [ obj ] = DCSeismicAnalysisBR()
            obj@DCBlakeRidge()
            
            obj.graphTop = obj.liuTop - obj.seafloorDepth; % mbsf
            obj.graphBottom = obj.liuBottom - obj.seafloorDepth; % mbsf
            
            obj.depthArray = (obj.minDepth : 1 : obj.maxDepth)';
            
            
            
            obj.SaturationLF = obj.LoadPhaseBehaviorBlakeRidge();
            
            
            
            
        end
        
        
        function [ SaturationLF ] = LoadPhaseBehaviorBlakeRidge( obj )
            % Loads saturation data from hydrate and gas grabit array .mat files and
            % saves into struct Saturation with fields 'Hydrate' and 'Gas'
            s = struct('Hydrate', cell2mat(struct2cell(load('Blake Ridge Data\Phase Behavior\HydrateBehavior.mat'))), ...
                       'Gas', cell2mat(struct2cell(load('Blake Ridge Data\Phase Behavior\GasBehavior.mat'))));

            % Upper and lower depths of grabit depth of investigation
            depthGrid = linspace(obj.liuTop, obj.liuBottom, 101);
            methaneQuantityGrid = linspace(0, 40, 41);
            
            [X, Y] = meshgrid(methaneQuantityGrid, depthGrid);
            
            % Hydrate saturations
            HydrateGrid = griddata(s.Hydrate(:,1), s.Hydrate(:,2), s.Hydrate(:,3), X, Y, 'cubic');
            HydrateGrid(HydrateGrid < 0) = 0;
            HydrateGrid(isnan(HydrateGrid)) = 0;
            HydrateGrid(HydrateGrid < .0001 & HydrateGrid ~= 0) = 0;
            
            % Gas saturations
            GasGrid = griddata(s.Gas(:,1), s.Gas(:,2), s.Gas(:,3), X, Y, 'cubic');
            GasGrid(GasGrid < 0) = 0;
            GasGrid(isnan(GasGrid)) = 0;
            GasGrid(GasGrid < .0001 & GasGrid ~= 0) = 0;
            
            % Hard coded fixes
            GasGrid(49, 40) = 0;
            GasGrid(66:end, 5) = GasGrid(end, 5);
            GasGrid(66:end, 6) = GasGrid(end, 6);
            GasGrid(65:end, 7) = GasGrid(end, 7);
            GasGrid(65:end, 8) = GasGrid(end, 8);
            GasGrid(65:end, 9) = GasGrid(end, 9);
            GasGrid(67:end, 10) = GasGrid(end, 10);
            
            quantityArray = (6:1:40)';
            n = numel(quantityArray);
            
            
            ThreePhaseTop = zeros(1, n);
            ThreePhaseBottom = zeros(1, n);
            
            HydrateFull = zeros(numel(obj.depthArray), n);
            GasFull = zeros(numel(obj.depthArray), n);
            
            
            for i = 1:n
                iQuantity = quantityArray(i);
                
                ThreePhaseTop(i) = sum(~any(GasGrid(:, iQuantity + 1), 2)) + 1;
                ThreePhaseBottom(i) = sum(any(HydrateGrid(:, iQuantity + 1), 2));
                
                
                HydrateFull(obj.saturationTop:obj.graphTop - 1, i) = HydrateGrid(1, i);
                HydrateFull(obj.graphTop:obj.graphBottom, i) = HydrateGrid(:, i);
                
                GasFull(obj.graphTop:obj.graphBottom, i) = GasGrid(:, i);
                GasFull(obj.graphBottom + 1:obj.saturationBottom, i) = GasGrid(end, i);        
            end
            
            
            
            SaturationLF.Top3P = ThreePhaseTop;
            SaturationLF.Bottom3P = ThreePhaseBottom;
            
            SaturationLF.Hydrate = HydrateFull;
            SaturationLF.Gas = GasFull;
            
            
            
        end
        function [ Dickens ] = LoadDickensBlakeRidge( depthArray , SaturationLF )
            % Returns CH4 quantity vs depth using Dickens et al. 1997
            % Input array is in mbsf
            % Output array is in rounded 

            % From grabit
            Dickens = load('Blake Ridge Data\Phase Behavior\MethaneQuantity.mat');
            % Column 1 is CH4 quantity
            % Colume 2 is mbsf
            Dickens.MethaneQuantity(:,1) = Dickens.MethaneQuantity(:,1).*16.04; % conversion from mol CH4 to g CH4 (per dm^3)

            % X is CH4 quantity
            % Y is mbsf
            [X,Y] = meshgrid( linspace(0,40,41) , depthArray );

            Sh_interpolated = interp2( X , Y , SaturationLF.Hydrate , Dickens.MethaneQuantity(:,1) , Dickens.MethaneQuantity(:,2) );
            Sg_interpolated = interp2( X , Y , SaturationLF.Gas , Dickens.MethaneQuantity(:,1) , Dickens.MethaneQuantity(:,2) );

            Dickens.Hydrate = interp1( Dickens.MethaneQuantity(:,2) , Sh_interpolated , depthArray );
            Dickens.Gas = interp1( Dickens.MethaneQuantity(:,2) , Sg_interpolated , depthArray );

            Dickens.Hydrate(isnan(Dickens.Hydrate)) = 0;
            Dickens.Gas(isnan(Dickens.Gas)) = 0;

            % finds average methane quantity over upper and lower bound of Dickens
            % upper_bound = 130;
            % lower_bound = 590;
            % interp_methane_quantity = interp1( Dickens.MethaneQuantity(:,2) , Dickens.MethaneQuantity(:,1) , depth_array( depth_array <= lower_bound & depth_array >= upper_bound ) );
            % interp_methane_quantity(isnan(interp_methane_quantity)) = [];
            % average_methane_quantity = mean(interp_methane_quantity);
            
            % Methane quantity vs depth
            figure

            axis_1 = subplot(2,1,1);
            h = plot( Dickens.MethaneQuantity(:,1) , Dickens.MethaneQuantity(:,2) , 'k' );
            set(gca,'YDir','Reverse')
            axis([0 40 100 600])
            % grid on
            xlabel('Methane quantity (g/dm^3 of pore volume)')
            ylabel('Depth (mbsf)')
            title('a) Experimentally Measured Methane Quantity')
            set(h,'LineWidth',2.5);

            % Fluid saturations vs depth
            % figure
            % scatter( Sh_interpolated , Dickens.MethaneQuantity(:,2) )
            % hold on
            % scatter( Sg_interpolated , Dickens.MethaneQuantity(:,2) )
            % hold on
            axis_2 = subplot(2,1,2);
            h = plot( Dickens.Hydrate , depthArray );
            set(h,'LineWidth',2.5);
            set(h,'Color',[0 .5 0]);
            hold on
            h = plot( Dickens.Gas , depthArray , 'r--' );
            set(h,'LineWidth',2.5);
            set(gca,'YDir','Reverse')
            axis([0 0.2 100 600])
            % axis([0 0.2 420 520])
            xlabel('Fluid saturations')
            ylabel('Depth (mbsf)')
            title('b) Variable Methane Quantity Scenario')
            % grid on
            legend('Hydrate saturation' , 'Gas saturation' )


            axis_1.Position = [.15 .59 .78 .37];
            axis_2.Position = [.15 .08 .78 .37];


            figure_1 = gcf;
            figure_1.Position(3) = 320;
            set(findall(figure_1,'-property','FontSize'),'FontSize',8)
            set(findall(figure_1,'-property','FontName'),'FontName','Arial')
        end
        
        
        
        
        
        
        
        
        
    end
    methods (Static)
        
    end
end

