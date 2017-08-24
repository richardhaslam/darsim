%  Basis functions updater F-AMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr & Matteo Cusini
%TU Delft
%Created: 07 August 2017
%Last modified: 14 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater_FAMS < bf_updater_ms
    properties
        Amedia
    end
    methods
        function ConstructPressureSystem(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections, Ntot)
            % Reservoir
            Km = ProductionSystem.Reservoir.K;
            S = ProductionSystem.CreateGlobalVariables(FineGrid, FluidModel.NofPhases, 'S_');
            Mob = FluidModel.ComputePhaseMobilities(S(:,1));
            Start = 1;
            End = FineGrid(1).N;
            obj.Amedia{1} = obj.MediumPressureSystem(FineGrid(1), Km, Mob(Start:End,:));
            obj.A = obj.Amedia{1};
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
        function MsP = MsProlongation(obj, FineGrid, CoarseGrid, Dimensions)
            % Prolongation operator for fractured reservoir
            MsP = [];
            Dimensions = Dimensions * ones(length(FineGrid), 1);
            Dimensions(2:end) = Dimensions(2:end) - 1;
            for i=1:length(FineGrid)
                cf = CoarseGrid(i).CoarseFactor ./ FineGrid(i).CoarseFactor;
                % Permutation Matrix
                [G, Ni, Nf, Ne, Nv] = obj.PermutationMatrix(FineGrid(i), CoarseGrid(i), cf);
                % Reorder A based on dual coarse grid partition 
                tildeA = G * obj.Amedia{i} * G';
                P = obj.CartesianMsP(tildeA, Ni, Nf, Ne, Nv, Dimensions(i));
                obj.Amedia{i} = P' * obj.Amedia{i} * P;
                MsP = blkdiag(MsP, G'*P);
            end          
        end
        function MsP = CartesianMsP(obj, tildeA, Ni, Nf, Ne, Nv, Dimensions)
            switch(Dimensions)
                case (1)
                    % 2. Define edge-edge (ee) block
                    Mee = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+1:Ni+Nf+Ne);
                    % 3. Define edge-node (en) block
                    Mev = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+Ne+1:Ni+Nf+Ne+Nv);
                    % 4. Compute inverse of (ii), (ff) and (ee) blocks
                    % 1D
                    Edges = slvblk(Mee, Mev);
                    
                    MsP = [-Edges;...
                        speye(Nv,Nv)];
                case (2)
                    % 2. Define face-face block
                    Mff = tildeA(Ni+1:Ni+Nf, Ni+1:Ni+Nf);
                    % 3. Define edge-face (fe) block
                    Mfe = tildeA(Ni+1:Ni+Nf, Ni+Nf+1:Ni+Nf+Ne);
                    % 4. Define edge-edge (ee) block
                    Mee = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+1:Ni+Nf+Ne) + diag(sum(Mfe,1));
                    % 5. Define edge-node (en) block
                    Mev = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+Ne+1:Ni+Nf+Ne+Nv);
                    % 6. Compute inverse of (ii), (ff) and (ee) blocks
                    % 2D
                    % when 2D there are no interiors.
                    Edges = slvblk(Mee, Mev);
                    Faces = slvblk(Mff, Mfe * Edges);
                    
                    MsP = [Faces;...
                        -Edges;...
                        speye(Nv,Nv)];
                case(3)
                    % 2. Define interior-interior (ii) block
                    Mii = tildeA(1:Ni, 1:Ni);
                    % 3. define interior-face (if) block
                    Mif = tildeA(1:Ni, Ni+1:Ni+Nf);
                    % 4. Define face-face block
                    Mff = tildeA(Ni+1:Ni+Nf, Ni+1:Ni+Nf) + diag(sum(Mif,1));
                    % Define edge-face (fe) block
                    Mfe = tildeA(Ni+1:Ni+Nf, Ni+Nf+1:Ni+Nf+Ne);
                    % 5. Define edge-edge (ee) block
                    Mee = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+1:Ni+Nf+Ne) + diag(sum(Mfe,1));
                    % 6. Define edge-node (en) block
                    Men = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+Ne+1:Ni+Nf+Ne+Nv);
                    % 3D
                    % 7. Compute inverse of (ii), (ff) and (ee) blocks
                    Edges = slvblk(Mee, Men);
                    Faces = slvblk(Mff, Mfe * Edges);
                    Interiors = slvblk(Mii, Mif * Faces);
                    
                    MsP = [-Interiors;...
                        Faces;...
                        -Edges;...
                        speye(Nv,Nv)];
            end
        end
    end
end