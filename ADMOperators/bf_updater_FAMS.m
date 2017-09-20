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
        BFtype = 1;
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
            % Prolongation operator for fractured reservoir (de-coupled)
            switch(obj.BFtype)
                case(1)
                    MsP = [];
                    Dimensions = Dimensions * ones(length(FineGrid), 1);
                    Dimensions(2:end) = Dimensions(2:end) - 1;
                    for i=1:length(FineGrid)
                        cf = CoarseGrid(i).CoarseFactor ./ FineGrid(i).CoarseFactor;
                        % Permutation Matrix
                        [G, Ni, Nf, Ne, Nv] = obj.PermutationMatrix(FineGrid(i), CoarseGrid(i), cf);
                        % Reorder A based on dual coarse grid partition
                        tildeA = G * obj.Amedia{i} * G';
                        P = obj.ComputeMsP(tildeA, Ni, Nf, Ne, Nv, Dimensions(i));
                        %obj.Amedia{i} = P' * obj.Amedia{i} * P;
                        MsP = blkdiag(MsP, G'*P);
                    end
                case(2)
                    cf = vertcat(CoarseGrid(1:end).CoarseFactor) ./ vertcat(FineGrid(1:end).CoarseFactor);
                    [G, Ni, Nf, Ne, Nv] = obj.FullyCoupledPermutationMatrix(FineGrid, CoarseGrid, cf);
                    tildeA = G * obj.A * G';
                    MsP = G'*obj.ComputeMsP(tildeA, Ni, Nf, Ne, Nv, Dimensions(1));
            end
        end
        function [P, nii, nff, nee, nvv] = FullyCoupledPermutationMatrix(obj, FineGrid, CoarseGrid, cf)
            %% 0. Define local variables 
            nxf = vertcat(FineGrid(1:end).Nx);
            nyf = vertcat(FineGrid(1:end).Ny);
            nzf = vertcat(FineGrid(1:end).Nz);
            nf  = vertcat(FineGrid(1:end).N );
            nxcf = cf(:,1);
            nycf = cf(:,2);
            nzcf = cf(:,3);
            
            nxc = vertcat(CoarseGrid(1:end).Nx);
            nyc = vertcat(CoarseGrid(1:end).Ny);
            nzc = vertcat(CoarseGrid(1:end).Nz);
            
            %% 1. Define number of cells of each block
            nii = ( (nxcf(1)-1)*(nycf(1)-1)*(nzcf(1)-1) ) * ( nxc(1) * nyc(1) * nzc(1) ); % interiors (only in matrix, no interiors in fractures)
            
            nff = ( (nxcf(1)-1)*(nycf(1)-1) + (nycf(1)-1)*(nzcf(1)-1) + (nxcf(1)-1)*(nzcf(1)-1) ) * nxc(1) * nyc(1) * nzc(1); % faces in matrix
            for f = 2:length(FineGrid)
                nff = nff + ( (nxcf(f)-1)*(nycf(f)-1) + (nycf(f)-1)*(nzcf(f)-1) + (nxcf(f)-1)*(nzcf(f)-1) ) * nxc(f) * nyc(f) * nzc(f); % faces in fractures
            end
            
            nee = ( (nxcf(1)-1) + (nycf(1)-1) + (nzcf(1)-1) ) * nxc(1) * nyc(1) * nzc(1); % edges in matrix
            for f = 2:length(FineGrid)
                nee = nee + ( (nxcf(f)-1) + (nycf(f)-1) + (nzcf(f)-1) ) * nxc(f) * nyc(f) * nzc(f); % edges in fractures
            end
            
            nvv = nxc(1) * nyc(1) * nzc(1); % verteces in matrix
            for f = 2:length(FineGrid)
                nvv = nvv + nxc(f) * nyc(f) * nzc(f); % verteces in fractures
            end
            
            P = sparse(sum(nf));
            
            %% 2. Construct matrix and vector
            % Coarsening factors are assumed to be odd numbers. 
            iii0 = 0;
            iff0 = nii;
            iee0 = nii + nff;
            ivv0 = nii + nff + nee;
            
            nxd = nxc + 1;
            nyd = nyc + 1;
            nzd = nzc + 1;
            icen = ceil(nxcf/2);
            jcen = ceil(nycf/2);
            kcen = ceil(nzcf/2);
            
            % Now fill in entries of permutation matrix 
            for m = 1:length(FineGrid)
                for k = 1:nzd(m)
                    for j = 1:nyd(m)
                        for i = 1:nxd(m)
                            for kz= 1:nzcf(m)
                                for jy= 1:nycf(m)
                                    for ix = 1:nxcf(m)
                                        % calculate the cell index in the original matrix
                                        ii = (i-1) * nxcf(m) - icen(m) + ix;
                                        jj = (j-1) * nycf(m) - jcen(m) + jy;
                                        kk = (k-1) * nzcf(m) - kcen(m) + kz;
                                        if ( ii >= 1 && ii <= nxf(m) && jj >= 1 && jj <= nyf(m) && kk >=1 && kk<=nzf(m))
                                            % There are dual domains that are
                                            % not withing the domain
                                            ijk = sum(nf(1:m-1)) + ii + (jj-1)*nxf(m) + (kk-1)*nxf(m)*nyf(m); % global index
                                            if (ix  ~= 1 && jy ~= 1 && kz ~=1)
                                                % internal point
                                                iii0 = iii0 + 1;
                                                P(iii0, ijk) = 1;
                                            elseif  (ix == 1 && jy == 1 && kz==1)
                                                % node points
                                                ivv0 = ivv0 + 1;
                                                P(ivv0, ijk) = 1;
                                            elseif (ix == 1 && jy == 1 && (kz ~=1 || kz~=nzcf(m)))
                                                % edge xy points
                                                iee0 = iee0 + 1;
                                                P(iee0, ijk) = 1;
                                            elseif (ix == 1 && kz == 1 && (jy ~=1 || jy~=nycf(m)))
                                                % edge xz points
                                                iee0 = iee0 + 1;
                                                P(iee0, ijk) = 1;
                                            elseif (jy == 1 && kz == 1 && (ix ~=1 || ix~=nxcf(m)))
                                                % edge yz points
                                                iee0 = iee0 + 1;
                                                P(iee0, ijk) = 1;
                                            else
                                                % face points
                                                iff0 = iff0 + 1;
                                                P(iff0, ijk) = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
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