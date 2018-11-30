%  Basis functions updater F-AMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr & Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater_FAMS < bf_updater_ms
    properties
        Amedia
        BFtype
    end
    methods
        function ConstructPressureSystem(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections)
            % Reservoir
            Km = ProductionSystem.Reservoir.K;
            S = ProductionSystem.CreateGlobalVariables(FineGrid, FluidModel.NofPhases, 'S_');
            Mob = FluidModel.ComputePhaseMobilities(S(:,1));
            Start = 1;
            End = FineGrid(1).N;
            obj.Amedia{1} = obj.MediumPressureSystem(FineGrid(1), Km, Mob(Start:End,:));
            obj.A = obj.Amedia{1};
            obj.AddWellsToPressureMatrix(ProductionSystem.Wells, Km, Mob(Start:End,:), FineGrid(1).N);
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Start = End + 1;
                End = Start + FineGrid(1+f).N - 1;
                Kf = ProductionSystem.FracturesNetwork.Fractures(f).K;
                obj.Amedia{1+f} = obj.MediumPressureSystem(FineGrid(1+f), Kf, Mob(Start:End,:));
                obj.A = blkdiag(obj.A, obj.Amedia{1+f});
            end
            % Non-Neighboring Connections
            Mobt = sum(Mob,2);
            for c = 1:length(CrossConnections)
                T_Geo = CrossConnections(c).T_Geo;
                i = c + FineGrid(1).N;
                j = CrossConnections(c).Cells;
                obj.A(i,j) = obj.A(i,j) - T_Geo' .* Mobt(j)';
                obj.A(i,i) = obj.A(i,i) + sum(T_Geo.* Mobt(j));
                obj.A(j,i) = obj.A(j,i) - T_Geo.* Mobt(i);
                obj.A(sub2ind(size(obj.A), j, j)) = obj.A(sub2ind(size(obj.A), j, j)) + T_Geo.* Mobt(i);
            end 
        end
        function [MsP, MsC] = MsProlongation(obj, FineGrid, CoarseGrid, Dimensions)
            % Prolongation operator for fractured reservoir (de-coupled)
            switch(obj.BFtype)
                case('DECOUPLED')
                    MsP = []; % Prolongation operator
                    MsC = []; % Correction functions operator
                    Dimensions = Dimensions * ones(length(FineGrid), 1);
                    Dimensions(2:end) = Dimensions(2:end) - 1;
                    for i=1:length(FineGrid)
                        cf = CoarseGrid(i).CoarseFactor ./ FineGrid(i).CoarseFactor;
                        % Permutation Matrix
                        [G, Ni, Nf, Ne, Nv] = obj.PermutationMatrix(FineGrid(i), CoarseGrid(i), cf);
                        % Reorder A based on dual coarse grid partition
                        tildeA = G * obj.Amedia{i} * G';
                        [P, C]= obj.ComputeMsP(tildeA, Ni, Nf, Ne, Nv, Dimensions(i));
                        %obj.Amedia{i} = P' * obj.Amedia{i} * P;
                        MsP = blkdiag(MsP, G'*P);
                        if i==1 && ~isempty(MsC)
                            MsC = blkdiag(MsC, G'*C*G);
                        end
                    end
                case('COUPLED')
                    MsP = []; % Prolongation operator
                    MsC = []; % Correction functions operator
                    cf = vertcat(CoarseGrid(1:end).CoarseFactor) ./ vertcat(FineGrid(1:end).CoarseFactor);
                    [G, Ni, Nf, Ne, Nv] = obj.FullyCoupledPermutationMatrix(FineGrid, CoarseGrid, cf);
                    tildeA = G * obj.A * G';
                    [MsP, MsC] = obj.ComputeMsP(tildeA, Ni, Nf, Ne, Nv, Dimensions(1));
                    MsP = G' * MsP;
                    if ~isempty(MsC)
                        MsC = G' * MsC * G;
                    end
            end
        end
        function UpdatePressureMatrix(obj, P, Grid)
            Start_f = 1;
            Start_c = 1;
            for m=1:length(obj.Amedia)
                [Nf, ~] = size(obj.Amedia{m});
                End_f = Start_f + Nf - 1;
                End_c = Start_c + Grid(m).N - 1 ;
                obj.Amedia{m} = P(Start_f:End_f, Start_c:End_c)' * obj.Amedia{m} * P(Start_f:End_f, Start_c:End_c);
                obj.Amedia{m} = obj.TransformIntoTPFA(obj.Amedia{m}, Grid(m));
                Start_f = End_f + 1;
                Start_c = End_c + 1;
            end
            obj.A = P' * obj.A * P;
        end
    end
end