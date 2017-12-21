%  Prolongation builder MS for hyperbolic variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 21 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_MSHyperbolic < prolongation_builder
    properties
        key = 'S_1';
        Pdelta
        Pdeltac
        nt;
    end
    methods
        function obj = prolongation_builder_MSHyperbolic(n)
            obj@prolongation_builder(n)
            obj.P = cell(1, n);
            obj.Pdeltac = cell(1, n);
            obj.nt = 1;
        end
        function BuildStaticOperators(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections, maxLevel, CoarseGrid)
            % Build Restriction and Prolongation operators for static grids
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(:, 1));
            obj.P{1} = obj.R{1}';
            for x = 2:maxLevel
                R = obj.MsRestriction(CoarseGrid(:, x-1), CoarseGrid(:,x));
                obj.R{x} = R*obj.R{x-1};
                obj.P{x} = obj.R{x}';
             end
        end
        function UpdateProlongationOperator(obj, FineGrid, CoarseGrid, ProductionSystem)
            %% Update the basis functions level by level (Update e = dSf/dSc)
            sf = ProductionSystem.Reservoir.State.Properties(obj.key).Value;
            sf0   =  ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
            for i=1:length(obj.R)
                obj.UpdateBasisFunctions(sf, sf0, i, CoarseGrid(i));
            end
        end
        function UpdateBasisFunctions(obj, sf, sf0, l, CoarseGrid)
            %% Update sat prolongation for level l
            % 1. Compute dSc from Sc^n = 1/V * sum(vSf^n)
            Sc  = obj.R{l}' * ((obj.R{l} * sf)  ./ sum(obj.R{l}, 2));
            Sc0 = obj.R{l}' * ((obj.R{l} * sf0) ./ sum(obj.R{l}, 2));
            deltac = Sc - Sc0;
            deltaf = sf - sf0;
            epsilon = deltaf ./ deltac;
            epsilon(isnan(epsilon)) = 1;
            epsilon(epsilon == Inf) = 1;
            epsilon(abs(deltac) < 1e-6) = 1;
            epsilon(abs(deltaf) < 1e-6) = 1;
%             if l==1
%                 cells = [10, 11, 12];
%                 Colors = {'ro', 'bo','go', 'blacko'};
%                 Names = {'cell 10 adm','cell 11 adm','cell 12 adm', 'dS coarse adm'};
%                 Line = {'r+', 'b+','g+', 'k+'};
% %                 Colors = {'r*', 'b*','g*', 'k*'};
% %                 Names = {'cell 900 fs','cell 1000 fs','cell 1100 fs', 'dS coarse fs'};
% %                 Line = {'r.', 'b.','g.', 'k.'};
%                 if obj.nt == 1
%                     figure(1)
%                     hold on
%                     for i=1:length(cells)
%                         plot(obj.nt, deltaf(cells(i)), Colors{i}, 'DisplayName', Names{i});
%                     end
%                     plot(obj.nt, deltac(cells(i)), Colors{i+1}, 'DisplayName', Names{i+1});
%                     xlabel('time-step');
%                     ylabel('deltaS');
%                     axis([0 210 0 5e-3]);
%                     legend('show');
%                     figure(2)
%                     hold on
%                     for i=1:length(cells)
%                         yyaxis left
%                         plot(obj.nt, epsilon(cells(i)), Colors{i}, 'DisplayName', Names{i})
%                     end
%                     xlabel('time-step');
%                     ylabel('epsilon');
%                     axis([0 210 0 2]);
%                     legend('show');
%                     for i=1:length(cells)
%                         yyaxis right
%                         plot(obj.nt, sf(cells(i)), Line{i}, 'DisplayName', Names{i});
%                     end
%                     plot(obj.nt, Sc(cells(i)), Line{i+1}, 'DisplayName', Names{i+1});
%                     xlabel('time-step');
%                     ylabel('Saturation');
%                     axis([0 210 0 1]);
%                     legend('show');
%                 else
%                     figure(1)
%                     hold on
%                     for i=1:length(cells)
%                         plot(obj.nt, deltaf(cells(i)), Colors{i}, 'HandleVisibility','off');
%                     end
%                     plot(obj.nt, deltac(cells(i)), Colors{i+1}, 'HandleVisibility','off')
%                     figure(2)
%                     hold on
%                     for i=1:length(cells)
%                         yyaxis left
%                         plot(obj.nt, epsilon(cells(i)), Colors{i}, 'HandleVisibility','off');
%                         yyaxis right
%                         plot(obj.nt, sf(cells(i)), Line{i}, 'HandleVisibility','off');
%                     end       
%                     plot(obj.nt, Sc(cells(i)), Line{i+1}, 'HandleVisibility','off'); 
%                 end                
%                 drawnow
%                 obj.nt = obj.nt+1;
%             end
            % 2. Update the prolongation operator
            obj.Pdelta  = deltaf;
            obj.Pdeltac{l} = deltac;
            for c=1:CoarseGrid.N
                fsI = CoarseGrid.GrandChildren(c,:);
                obj.P{l}(fsI, c) = epsilon(fsI);
                %obj.P{l}(fsI, c) = epsilon(fsI)/max(epsilon(fsI));
            end
        end
        function ADMProl = ADMProlongation(obj, ADMGrid, GlobalGrids, ADMRest)
            ADMProl = ADMRest';
            % Coarse levels cells
            for c = ADMGrid.N(1) + 1:ADMGrid.Ntot
                indexes = ADMGrid.GrandChildren{c};
                ADMProl(indexes, c) = obj.P{ADMGrid.level(c)}(indexes, ADMGrid.CellIndex(c));
            end
        end
        function AverageMassOnCoarseBlocks(obj, Formulation, ProductionSystem, FineGrid, FluidModel, ADMRest)
            % virtual call
        end
    end
end