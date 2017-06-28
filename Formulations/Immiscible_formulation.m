% Immiscible Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 7 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_formulation < formulation
    properties
        % Sequential run variables
        Mobt
        Utot
        Qwells
        f
        df
        V
    end
    methods
        function obj = Immiscible_formulation()
            obj@formulation();
            obj.Tph = cell(2,1);
            obj.Gph = cell(2,1);
        end
        function x = GetPrimaryUnknowns(obj, ProductionSystem, DiscretizationModel)
             Nt = DiscretizationModel.N;
             Nm = DiscretizationModel.ReservoirGrid.N;
             x = zeros(obj.NofPhases * Nt, 1);
             %% Reservoir
             Start = 1;
             End = Nm;
             x(Start:End) = ProductionSystem.Reservoir.State.Properties('P_2').Value;
             for i=1:obj.NofPhases-1
                 Start = End + 1;
                 End = Start + Nm - 1;
                 x(Start:End) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(i)]).Value;
             end
             %% Fractures
             if ProductionSystem.FracturesNetwork.Active
                 Nf = DiscretizationModel.FracturesGrid.N;
                 for f=1:ProductionSystem.FracturesNetwork.NumOfFrac
                     Start = End + 1;
                     End = Start + Nf(f) - 1;
                     x(Start:End) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_2').Value;
                     for i=1:obj.NofPhases-1
                         Start = End + 1;
                         End = Start + Nf(f) - 1;
                         x(Start:End) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['S_', num2str(i)]).Value;
                     end
                 end
             end
        end
        function x = GetPrimaryPressure(obj, ProductionSystem, DiscretizationModel)
            Nt = DiscretizationModel.N;
            Nm = DiscretizationModel.ReservoirGrid.N;
            %Nf = DiscretizationModel.FracturesGrid.N;
            x = zeros(Nt, 1);
            if obj.NofPhases > 1
                x(1:Nm) = ProductionSystem.Reservoir.State.Properties('P_2').Value;
            else
                x(1:Nm) = ProductionSystem.Reservoir.State.Properties('P_1').Value;
            end
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            %% 1. Reservoir Properteis and Derivatives
            obj.drhodp = FluidModel.DrhoDp(ProductionSystem.Reservoir.State);
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.dPc = FluidModel.DPcDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            %% 2. Fractures Properteis and Derivatives
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                obj.drhodp = [obj.drhodp; FluidModel.DrhoDp(ProductionSystem.FracturesNetwork.Fractures(f).State) ];
                obj.Mob = [obj.Mob; FluidModel.ComputePhaseMobilities(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('S_1').Value)];
                obj.dMob = [obj.dMob; FluidModel.DMobDS(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('S_1').Value)];
                obj.dPc = [obj.dPc; FluidModel.DPcDS(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('S_1').Value)];
            end
        end
        %% Methods for FIM Coupling
        function Residual = BuildMediumResidual(obj, Medium, Grid, dt, State0, Index, qw, qf, f)
            % Create local variables
            N = Grid.N;
            pv = Medium.Por*Grid.Volume;
            s_old = zeros(N, obj.NofPhases);
            rho_old = zeros(N, obj.NofPhases);
            P = zeros(N, obj.NofPhases);
            s = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            
            % Copy values in local variables
            for i=1:obj.NofPhases
                s_old(:,i) = State0.Properties(['S_', num2str(i)]).Value(Index.Start:Index.End);
                rho_old(:,i) = State0.Properties(['rho_', num2str(i)]).Value(Index.Start:Index.End);
                P(:, i) = Medium.State.Properties(['P_', num2str(i)]).Value;
                s(:, i) = Medium.State.Properties(['S_', num2str(i)]).Value;
                rho(:, i) = Medium.State.Properties(['rho_', num2str(i)]).Value;
            end
            depth = Grid.Depth;
            
            % Accumulation Term
            AS = speye(N)*pv/dt;

            % RESIDUAL
            Residual = zeros(obj.NofPhases*N,1);
            for i=1:obj.NofPhases
                Residual((i-1)*N+1:i*N)  = AS*(rho(:,i) .* s(:,i) - rho_old(:,i) .* s_old(:,i))...
                                           + obj.Tph{i, 1+f} * P(:, i)...
                                           - obj.Gph{i, 1+f} * depth...
                                           - qw(Index.Start:Index.End, i)...
                                           - qf(Index.Start:Index.End, i);
            end
        end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute vector of qs
            qw = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            qf = zeros(DiscretizationModel.N, obj.NofPhases);
            if ProductionSystem.FracturesNetwork.Active
                qf = ComputeSourceTerms_frac_mat(obj, ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures, DiscretizationModel);
            end
            % Transmissibility of reservoir
            for i=1:obj.NofPhases
                [obj.Tph{i, 1}, obj.Gph{i, 1}] = ...
                    obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, obj.UpWind{i, 1}, obj.Mob(1:DiscretizationModel.ReservoirGrid.N, i), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(i)]).Value, obj.GravityModel.RhoInt{i, 1});
            end
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            Nm = DiscretizationModel.ReservoirGrid.N;
            % BuildResidual for Reservoir
            Index.Start = 1;
            Index.End = Nm;
            Residual = BuildMediumResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, qw, qf, 0);
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                Index.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, DiscretizationModel.FracturesGrid.N(f));
                % Transmissibility of fractures cells
                for i=1:obj.NofPhases
                    [obj.Tph{i, 1+f}, obj.Gph{i, 1+f}] = ...
                        obj.TransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), obj.UpWind{i, 1+f}, obj.Mob(Index.Start:Index.End, i),...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(i)]).Value, obj.GravityModel.RhoInt{i, 1+f});
                end
                % BuildResidual for Fractures
                Residual_frac_f = BuildMediumResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, qw, qf, f);
                Residual = [Residual; Residual_frac_f];
            end
        end
        function JacInd   = DefineJacobianIndex(obj, DiscretizationModel)
            JacInd = struct;
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nf = DiscretizationModel.FracturesGrid.N;
            for If1_Local = 1 : length(DiscretizationModel.CrossConnections)
                Cells = DiscretizationModel.CrossConnections(If1_Local).Cells;
                
                If1_Global = If1_Local + Nm;
                Index_frac1_Local = DiscretizationModel.Index_Global_to_Local(If1_Global);
                indices_m = Cells( Cells <= Nm );
                for i = 1:obj.NofPhases
                    % frac-matrix
                    % Pressure Blocks
                    % Adding the Jp_frac_mat of phase i to Jacobian
                    % adding to Jp_FM
                    row_ind_Jp{1,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                    col_ind_Jp{1,i} = indices_m;
                    
                    % adding to Jp_MF
                    col_ind_Jp{2,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + Index_frac1_Local.g;
                    row_ind_Jp{2,i} = (i-1)*Nm + indices_m;
                    
                    % Saturation Blocks 
                    % Adding the Js_frac_mat of phase i to matrix columns
                    % adding to Js_FM
                    row_ind_Js{1,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                    col_ind_Js{1,i} = indices_m + Nm;
            
                    % Adding the Js_frac_mat of phase i to fractures columns
                    % adding to Js_MF
                    row_ind_Js{2,i} = (i-1)*Nm + indices_m;
                    col_ind_Js{2,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + Nf(Index_frac1_Local.f) + Index_frac1_Local.g;
                    
                    % adding to Js_FF
                    row_ind_Js{3,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                    col_ind_Js{3,i} = [];
                    
                    
                    % frac-frac
                    indices_f = DiscretizationModel.CrossConnections(If1_Local).Cells( DiscretizationModel.CrossConnections(If1_Local).Cells > Nm );
                    If2_Global = [];
                    If2_Local = [];
                    Index_frac2_Local = [];
                    T_Geo_Harmo = [];
                    if ~isempty(indices_f)
                        for n = 1:length(indices_f)
                            If2_Global{n} = indices_f(n); % Global indices of the other fractures' cells if any
                            If2_Local{n} = If2_Global - Nm;
                            Index_frac2_Local{n} = DiscretizationModel.Index_Global_to_Local(If2_Global);
                            T_Geo_Half1 = DiscretizationModel.CrossConnections(If1_Local).T_Geo( length(indices_m)+n );
                            T_Geo_Half2 = DiscretizationModel.CrossConnections(If2_Local).T_Geo( DiscretizationModel.CrossConnections(If2_Local).Cells==If1_Global );
                            T_Geo_Harmo{n} = (T_Geo_Half1 * T_Geo_Half2) / (T_Geo_Half1 + T_Geo_Half2);
                            
                            % Pressure Blocks
                            % Adding the Jp_frac_frac of phase i to Jacobian
                            % adding to Jp_F1F2
                            row_ind_Jp{3,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                            col_ind_Jp{3,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) +                                 + Index_frac2_Local.g;
                            
                            % adding to Jp_F1F1
                            row_ind_Jp{4,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                            col_ind_Jp{4,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) +                                 + Index_frac1_Local.g;
                            
                            % adding to Jp_F2F1
                            row_ind_Jp{5,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) + (i-1)*(Nf(Index_frac2_Local.f)) + Index_frac2_Local.g;
                            col_ind_Jp{5,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) +                                 + Index_frac1_Local.g;
                            
                            % adding to Jp_F2F2
                            row_ind_Jp{6,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) + (i-1)*(Nf(Index_frac2_Local.f)) + Index_frac2_Local.g;
                            col_ind_Jp{6,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) +                                 + Index_frac2_Local.g;
                                       
                            % Saturation Blocks
                            % Adding the Js_frac2_frac1 of phase i to frac2 columns
                            % adding to Js_F1F2
                            row_ind_Js{4,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                            col_ind_Js{4,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) +    1 *(Nf(Index_frac2_Local.f)) + Index_frac2_Local.g;
                            
                            % adding to Js_F2F2
                            row_ind_Js{5,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) + (i-1)*(Nf(Index_frac2_Local.f)) + Index_frac2_Local.g;
                            col_ind_Js{5,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) +     1*(Nf(Index_frac2_Local.f)) + Index_frac2_Local.g;
                            
                            % Adding the Js_frac1_frac2 of phase i to frac1 columns
                            % adding to Js_F2F1
                            row_ind_Js{6,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) + (i-1)*(Nf(Index_frac2_Local.f)) + Index_frac2_Local.g;
                            col_ind_Js{6,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) +    1 *(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                            
                            % adding to Js_F1F1
                            row_ind_Js{7,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                            col_ind_Js{7,i} = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) +    1 *(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                            
                        end
                    end
                end
                JacInd(If1_Local).If1_Global = If1_Global;
                JacInd(If1_Local).Index_frac1_Local = Index_frac1_Local;
                JacInd(If1_Local).indices_m = indices_m;
                JacInd(If1_Local).row_ind_Jp = row_ind_Jp;
                JacInd(If1_Local).col_ind_Jp = col_ind_Jp;
                JacInd(If1_Local).row_ind_Js = row_ind_Js;
                JacInd(If1_Local).col_ind_Js = col_ind_Js;
                JacInd(If1_Local).indices_f = indices_f;
                JacInd(If1_Local).If2_Global = If2_Global;
                JacInd(If1_Local).If2_Local = If2_Local;
                JacInd(If1_Local).Index_frac2_Local = Index_frac2_Local;
                JacInd(If1_Local).T_Geo_Harmo = T_Geo_Harmo;
            end
        end
        function Jacobian = BuildMediumJacobian(obj, Medium, Wells, Grid, dt, Index, f)
            % Create local variables
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            pv = Grid.Volume*Medium.Por;
            rho = zeros(N, obj.NofPhases);
            s = zeros(N, obj.NofPhases);
            for i=1:obj.NofPhases
                rho(:,i) = Medium.State.Properties(['rho_', num2str(i)]).Value;
                s(:,i) = Medium.State.Properties(['S_', num2str(i)]).Value;
            end
            
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            Jp = cell(2,1);
            JS = cell(2,1);
            
            for i=1:obj.NofPhases
                % 1.a Pressure Block
                Jp{i} = obj.Tph{i,1+f};

                % 1.b: compressibility part
                dMupx = obj.UpWind{i,1+f}.x*(obj.Mob(Index.Start:Index.End, i) .* obj.drhodp(Index.Start:Index.End,i));
                dMupy = obj.UpWind{i,1+f}.y*(obj.Mob(Index.Start:Index.End, i) .* obj.drhodp(Index.Start:Index.End,i));
                dMupz = obj.UpWind{i,1+f}.z*(obj.Mob(Index.Start:Index.End, i) .* obj.drhodp(Index.Start:Index.End,i));
                
                vecX1 = min(reshape(obj.U{i,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{i,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{i,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{i,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{i,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{i,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz; 
                acc = pv/dt .* obj.drhodp(Index.Start:Index.End,i) .* s(:,i);
                
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                Jp{i} = Jp{i} + spdiags(DiagVecs, DiagIndx, N, N);
                
                % 2. Saturation Block
                dMupx = obj.UpWind{i,1+f}.x * (obj.dMob(Index.Start:Index.End, i) .* rho(:,i));
                dMupy = obj.UpWind{i,1+f}.y * (obj.dMob(Index.Start:Index.End, i) .* rho(:,i));
                dMupz = obj.UpWind{i,1+f}.z * (obj.dMob(Index.Start:Index.End, i) .* rho(:,i));
                % Construct JS block
                x1 = min(reshape(obj.U{i,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                x2 = max(reshape(obj.U{i,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                y1 = min(reshape(obj.U{i,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                y2 = max(reshape(obj.U{i,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                z1 = min(reshape(obj.U{i,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                z2 = max(reshape(obj.U{i,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
                v = (-1)^(i+1) * ones(N,1)*pv/dt .* rho(:,i);
                DiagVecs = [-z2, -y2, -x2, z2+y2+x2-z1-y1-x1+v, x1, y1, z1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                JS{i} = spdiags(DiagVecs,DiagIndx,N,N);
            end  
            
            %Add capillarity
            JS {1}= JS{1} - Jp{1} * spdiags(obj.dPc, 0, N, N);
            
            % Add Wells
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [Jp{1}, JS{1}, Jp{2}, JS{2}] = obj.AddWellsToJacobian(Jp{1}, JS{1}, Jp{2}, JS{2}, Medium.State, Wells, Medium.K(:,1));
            end
            % Full Jacobian: put the 4 blocks together
            Jacobian = [Jp{1}, JS{1}; Jp{2}, JS{2}];
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %    pm  sm  | pf   sf
            % m| Jp1 Js1 | Jp1 Js1 |
            % m| Jp2 Js2 | Jp2 Js2 |
            %  |---------|---------|
            % f| Jp1 Js1 | Jp1 Js1 |
            % f| Jp2 Js2 | Jp2 Js2 |
            %  |---------|---------|
            
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            Reservoir = ProductionSystem.Reservoir;
            Fractures = ProductionSystem.FracturesNetwork.Fractures;
            Wells = ProductionSystem.Wells;
            
            % Jacobian of the reservoir
            Index.Start = 1;
            Index.End = Nm;
%             start1 = tic;
            Jacobian = BuildMediumJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0);
            
            % Jacobian of the fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Nf = DiscretizationModel.FracturesGrid.N;
                Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                Index.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, DiscretizationModel.FracturesGrid.Grids(f).N);
                Jacobian_frac_f = BuildMediumJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f);
                Jacobian  = blkdiag( Jacobian , Jacobian_frac_f );
            end
%             disp(['Building Jacobian without cross conn took: ', num2str(toc(start1))]);
            
            % Effects of frac-matrix & frac-frac connections in Jacobian
%             start = tic;
            JacInd =  DiscretizationModel.JacInd;
            for If1_Local = 1 : length(DiscretizationModel.CrossConnections)
                DiscretizationModel.CrossConnections(If1_Local).CrossConnectionsUpWind(ProductionSystem, DiscretizationModel, If1_Local, obj);
                Cells = DiscretizationModel.CrossConnections(If1_Local).Cells;
                T_Geo = DiscretizationModel.CrossConnections(If1_Local).T_Geo;
                UpWind = DiscretizationModel.CrossConnections(If1_Local).UpWind;
                
                If1_Global = JacInd(If1_Local).If1_Global;
                Index_frac1_Local = JacInd(If1_Local).Index_frac1_Local;
                indices_m = JacInd(If1_Local).indices_m;
                indices_f = JacInd(If1_Local).indices_f;
                If2_Global = JacInd(If1_Local).If2_Global;
                If2_Local = JacInd(If1_Local).If2_Local;
                Index_frac2_Local = JacInd(If1_Local).Index_frac2_Local;
                T_Geo_Harmo = JacInd(If1_Local).T_Geo_Harmo;
                
                for i = 1:obj.NofPhases
                    % Loop over equations
                    
                    rho_f1 = Fractures(Index_frac1_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac1_Local.g);
                    Pf1 = Fractures(Index_frac1_Local.f).State.Properties(['P_',num2str(i)]).Value(Index_frac1_Local.g);

                    % frac-matrix
                    % Pressure Blocks
                    
                    rho_m = Reservoir.State.Properties(['rho_',num2str(i)]).Value(indices_m);
                    Pm = Reservoir.State.Properties(['P_',num2str(i)]).Value(indices_m);
                    Jp_frac_mat = - T_Geo(1:length(indices_m)) .* (...
                         UpWind(1:length(indices_m),i).*obj.Mob(indices_m,i) .*( rho_m  + obj.drhodp(indices_m,i) .*(Pm-Pf1) ) + ...
                        ~UpWind(1:length(indices_m),i).*obj.Mob(If1_Global,i).*( rho_f1 + obj.drhodp(If1_Global,i).*(Pm-Pf1) ) );
                    
                    row_ind_Jp = JacInd(If1_Local).row_ind_Jp;
                    col_ind_Jp = JacInd(If1_Local).col_ind_Jp;
                    row_ind_Js = JacInd(If1_Local).row_ind_Js;
                    col_ind_Js = JacInd(If1_Local).col_ind_Js;
                    % Adding the Jp_frac_mat of phase i to Jacobian
                    % adding to Jp_FM
                    Jacobian(row_ind_Jp{1,i}, col_ind_Jp{1,i}) = Jp_frac_mat';
                    % adding to Jp_FF
                    Jacobian(row_ind_Jp{1,i}, row_ind_Jp{1,i} - (i-1)*Nf(Index_frac1_Local.f)) = ...
                    Jacobian(row_ind_Jp{1,i}, row_ind_Jp{1,i} - (i-1)*Nf(Index_frac1_Local.f)) - sum(Jp_frac_mat);
                    % adding to Jp_MF
                    Jacobian(row_ind_Jp{2,i}, col_ind_Jp{2,i}) = Jp_frac_mat;
                    % adding to Jp_MM
                    Jacobian(sub2ind([Nt*obj.NofPhases-(i-1)*Nm , Nm], indices_m, indices_m) + (i-1).*Nm.*indices_m) = ...
                    Jacobian(sub2ind([Nt*obj.NofPhases-(i-1)*Nm , Nm], indices_m, indices_m) + (i-1).*Nm.*indices_m) - Jp_frac_mat;
  
                    % Saturation Blocks
                    Js_mat_frac = T_Geo(1:length(indices_m)) .*  UpWind(1:length(indices_m),i).*(Pf1-Pm).*( rho_m .* obj.dMob(indices_m ,i) );
                    Js_frac_mat = T_Geo(1:length(indices_m)) .* ~UpWind(1:length(indices_m),i).*(Pf1-Pm).*( rho_f1.* obj.dMob(If1_Global,i) );
                    % Adding the Js_frac_mat of phase i to matrix columns
                    % adding to Js_FM
                    Jacobian(row_ind_Js{1,i}, col_ind_Js{1,i}) = Js_mat_frac;
                    % adding to Js_MM
                    Jacobian(sub2ind([Nt*obj.NofPhases-(i-1)*Nm , Nm], indices_m, indices_m) + obj.NofPhases*Nm*Nt + (i-1).*Nm.*indices_m) = ...
                    Jacobian(sub2ind([Nt*obj.NofPhases-(i-1)*Nm , Nm], indices_m, indices_m) + obj.NofPhases*Nm*Nt + (i-1).*Nm.*indices_m) - Js_mat_frac;
                    % Adding the Js_frac_mat of phase i to fractures columns
                    % adding to Js_MF
                    Jacobian(row_ind_Js{2,i}, col_ind_Js{2,i}) = -Js_frac_mat;
                    % adding to Js_FF
                    Jacobian(row_ind_Js{3,i}, row_ind_Js{3,i} + (obj.NofPhases-i)*Nf(Index_frac1_Local.f)) = ...
                    Jacobian(row_ind_Js{3,i}, row_ind_Js{3,i} + (obj.NofPhases-i)*Nf(Index_frac1_Local.f)) + sum(Js_frac_mat);
                    
                    % frac-frac
                    if ~isempty(indices_f)
                        for n = 1:length(indices_f)
                            If2_Global{n} = indices_f(n); % Global indices of the other fractures' cells if any
                            If2_Local{n} = If2_Global - Nm;
                            Index_frac2_Local{n} = DiscretizationModel.Index_Global_to_Local(If2_Global);
                            rho_f2 = Fractures(Index_frac2_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac2_Local.g);
                            Pf2 = Fractures(Index_frac2_Local.f).State.Properties(['P_',num2str(i)]).Value(Index_frac2_Local.g);
                            T_Geo_Harmo{n} = (T_Geo_Half1 * T_Geo_Half2) / (T_Geo_Half1 + T_Geo_Half2);
                            
                            % Pressure Blocks
                            Jp_frac_frac = - T_Geo_Harmo * (...
                                UpWind(length(indices_m)+n,i) * obj.Mob(If1_Global,i) * ( rho_f2 + obj.drhodp(If2_Global,i)*(Pf2-Pf1) ) + ...
                               ~UpWind(length(indices_m)+n,i) * obj.Mob(If2_Global,i) * ( rho_f1 + obj.drhodp(If1_Global,i)*(Pf2-Pf1) ) );
                            % Adding the Jp_frac_frac of phase i to Jacobian
                            % adding to Jp_F1F2
                            Jacobian(row_ind_Jp{3,i}, col_ind_Jp{3,i}) = Jp_frac_frac;
                            % adding to Jp_F1F1
                            Jacobian(row_ind_Jp{4,i}, col_ind_Jp{4,i}) = Jacobian(row_ind_Jp, col_ind_Jp) - Jp_frac_frac;
                            % adding to Jp_F2F1
                            Jacobian(row_ind_Jp{5,i}, col_ind_Jp{5,i}) = Jp_frac_frac;
                            % adding to Jp_F2F2
                            Jacobian(row_ind_Jp{6,i}, col_ind_Jp{6,i}) = Jacobian(row_ind_Jp, col_ind_Jp) - Jp_frac_frac;
                           
                            % Saturation Blocks
                            Js_frac2_frac1 = T_Geo_Harmo* UpWind(length(indices_m)+n,i)*(Pf1-Pf2)*(rho_f2*obj.dMob(If2_Global,i));
                            Js_frac1_frac2 = T_Geo_Harmo*~UpWind(length(indices_m)+n,i)*(Pf1-Pf2)*(rho_f1*obj.dMob(If1_Global,i));

                            % Adding the Js_frac2_frac1 of phase i to frac2 columns
                            % adding to Js_F1F2
                            Jacobian(row_ind_Js{4,i}, col_ind_Js{4,i}) = Js_frac2_frac1;
                            % adding to Js_F2F2
                            Jacobian(row_ind_Js{5,i}, col_ind_Js{5,i}) = Jacobian(row_ind_Js, col_ind_Js) - Js_frac2_frac1;
                            % Adding the Js_frac1_frac2 of phase i to frac1 columns
                            % adding to Js_F2F1
                            Jacobian(row_ind_Js{6,i}, col_ind_Js{6,i}) = Js_frac1_frac2;
                            % adding to Js_F1F1
                            Jacobian(row_ind_Js{7,i}, col_ind_Js{7,i}) = Jacobian(row_ind_Js, col_ind_Js) - Js_frac1_frac2;
                        end
                    end
                end
%                 disp(['Time for 1 frac: ', num2str(toc(onefrac))]);
            end
%             Jacobian = sparse(rows, columns, values, size, size) + Jacobian;
%             disp(['Adding cross connection to the Jacobian took ',  num2str(toc(start))]);
        end
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                Nm = DiscretizationModel.ReservoirGrid.N;
                %% 1. Update matrix
                % Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(obj.NofPhases)]);
                Pm.update(delta(1:Nm));
                DeltaLast = zeros(Nm, 1);
                for ph = 1:obj.NofPhases-1
                    Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(ph)]);
                    DeltaS = delta(ph*Nm + 1:(ph+1)*Nm);
                    Sm.update(DeltaS);
                    % Remove values that are not physical
                    Sm.Value = max(Sm.Value, 0);
                    Sm.Value = min(Sm.Value, 1);
                    DeltaLast = DeltaLast + DeltaS;
                end
                Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(obj.NofPhases)]);
                Sm.update(-DeltaLast);
                % Remove values that are not physical
                Sm.Value = max(Sm.Value, 0);
                Sm.Value = min(Sm.Value, 1);
                % Update Phase Densities
                FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
                % Update total density
                FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State);
                % Update Pc
                FluidModel.ComputePc(ProductionSystem.Reservoir.State);

                %% 2. Update fractures pressure and densities
                if ProductionSystem.FracturesNetwork.Active
                    Nf = DiscretizationModel.FracturesGrid.N;
                    for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                        % Update Pressure
                        Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(obj.NofPhases)]);
                        Pf.update(delta( obj.NofPhases*(Nm+sum(Nf(1:f-1)))+1 : obj.NofPhases*(Nm+sum(Nf(1:f-1)))+Nf(f) ));
                        DeltaLast = zeros(Nf(f), 1);
                        for ph = 1:obj.NofPhases-1
                            Sf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['S_', num2str(ph)]);
                            Sf.update(delta( obj.NofPhases*(Nm+sum(Nf(1:f-1)))+ ph*Nf(f) + 1 : obj.NofPhases*(Nm+sum(Nf(1:f-1))) + (ph+1)*Nf(f)));
                            Sf.Value = max(Sf.Value, 0);
                            Sf.Value = min(Sf.Value, 1);
                            DeltaLast = DeltaLast + delta(obj.NofPhases*(Nm+sum(Nf(1:f-1)))+ ph*Nf(f) + 1 : obj.NofPhases*(Nm+sum(Nf(1:f-1))) + (ph+1)*Nf(f));
                        end
                        Sf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['S_', num2str(obj.NofPhases)]);
                        Sf.update(-DeltaLast);
                        Sf.Value = max(Sf.Value, 0);
                        Sf.Value = min(Sf.Value, 1);
                        % Update Phase Densities
                        FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                        % Update total density
                        FluidModel.ComputeTotalDensity(ProductionSystem.FracturesNetwork.Fractures(f).State);
                        % Update Pc
                        FluidModel.ComputePc(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    end
                end
            end
        end
        function [Tph, Gph] = TransmissibilityMatrix(obj, Grid, UpWind, Mob, rho, RhoInt)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            % Transmissibility matrix construction
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            
            % Apply upwind operator
            Mupx = UpWind.x*(Mob .* rho);
            Mupy = UpWind.y*(Mob .* rho);
            Mupz = UpWind.z*(Mob .* rho);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);
            Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1);
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1);
            y1 = reshape(Ty(:,1:Ny,:), N, 1);
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1);
            z1 = reshape(Tz(:,:,1:Nz), N, 1);
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1);
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            Tph = spdiags(DiagVecs,DiagIndx,N,N);
            
            % Gravity Matrix
            Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* RhoInt.x(2:Grid.Nx,:,:);
            Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* RhoInt.y(:,2:Grid.Ny,:);
            Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* RhoInt.z(:,:,2:Grid.Nz);
            
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1);
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1);
            y1 = reshape(Ty(:,1:Ny,:), N, 1);
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1);
            z1 = reshape(Tz(:,:,1:Nz), N, 1);
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1);
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            Gph = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function q = ComputeSourceTerms(obj, N, Wells)
            q = zeros(N, obj.NofPhases);    
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj(i).Cells;
                q(c, :) = Wells.Inj(i).QPhases(:,:);
            end
            %Producers
            for i=1:Wells.NofProd
                c = Wells.Prod(i).Cells;
                q(c, :) = Wells.Prod(i).QPhases(:,:);
            end
        end
        function [qf] = ComputeSourceTerms_frac_mat(obj, Reservoir, Fractures, DiscretizationModel)
            qfm = zeros(DiscretizationModel.ReservoirGrid.N, obj.NofPhases);
            qmf = zeros(sum(DiscretizationModel.FracturesGrid.N), obj.NofPhases);
            for If1_Local = 1 : length(DiscretizationModel.CrossConnections)
                indices_m = DiscretizationModel.CrossConnections(If1_Local).Cells( DiscretizationModel.CrossConnections(If1_Local).Cells <= DiscretizationModel.ReservoirGrid.N );
                If1_Global = DiscretizationModel.ReservoirGrid.N+If1_Local; % Global index of this fracture cells;
                Index_frac1_Local = DiscretizationModel.Index_Global_to_Local(If1_Global);
                for i = 1 : obj.NofPhases
                    % matrix-frac & frac-matrix
                    qmf(If1_Local, i) = qmf(If1_Local, i) + ...
                        sum(...
                        DiscretizationModel.CrossConnections(If1_Local).T_Geo(1:length(indices_m)) * obj.Mob(If1_Global, i) .* ...
                        Reservoir.State.Properties(['rho_', num2str(i)]).Value(indices_m) .*...
                        (Reservoir.State.Properties(['P_', num2str(i)]).Value(indices_m) - Fractures(Index_frac1_Local.f).State.Properties(['P_', num2str(i)]).Value( Index_frac1_Local.g))...
                        );
                    qfm(indices_m, i) = qfm(indices_m, i) + ...
                        DiscretizationModel.CrossConnections(If1_Local).T_Geo(1:length(indices_m)) .* obj.Mob(If1_Global, i) .* ...
                        Reservoir.State.Properties(['rho_', num2str(i)]).Value(indices_m) .*...
                        (Fractures(Index_frac1_Local.f).State.Properties(['P_', num2str(i)]).Value(Index_frac1_Local.g) - Reservoir.State.Properties(['P_', num2str(i)]).Value(indices_m));
                    
                    % frac1-frac2 & frac2-frac1
                    indices_f = DiscretizationModel.CrossConnections(If1_Local).Cells( DiscretizationModel.CrossConnections(If1_Local).Cells > DiscretizationModel.ReservoirGrid.N );
                    if ~isempty(indices_f)
                        for n = 1:length(indices_f) 
                            If2_Global = indices_f(n); % Global indices of the other fractures' cells if any
                            If2_Local = If2_Global - DiscretizationModel.ReservoirGrid.N;
                            T_Geo_Half1 = DiscretizationModel.CrossConnections(If1_Local      ).T_Geo( length(indices_m)+n );
                            T_Geo_Half2 = DiscretizationModel.CrossConnections(If2_Local).T_Geo( DiscretizationModel.CrossConnections(If2_Local).Cells==If1_Global );
                            T_Geo_Harmo = (T_Geo_Half1 * T_Geo_Half2) / (T_Geo_Half1 + T_Geo_Half2);
                            Index_frac2_Local = DiscretizationModel.Index_Global_to_Local(If2_Global);
                            
                            qmf(If1_Local, i) = qmf(If1_Local, i) + T_Geo_Harmo * obj.Mob(If1_Global, i) * ...
                                Fractures(Index_frac1_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac1_Local.g) * ...
                                (Fractures(Index_frac2_Local.f).State.Properties(['P_', num2str(i)]).Value(Index_frac2_Local.g) - ...
                                 Fractures(Index_frac1_Local.f).State.Properties(['P_', num2str(i)]).Value(Index_frac1_Local.g));
                             % As the whole loop is going through all the fractures cells, there is no need to write a 2nd line of code to ...
                             % add the qmf term for the frac2-frac1 connection.
                        end
                    end
                end
            end
            qf = [qfm;qmf];
        end
        function [Jwp, JwS, Jnwp, JnwS] = AddWellsToJacobian(obj, Jwp, JwS, Jnwp, JnwS, State, Wells, K)
            % Define Local handles
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            p = State.Properties('P_2').Value;
            rho = zeros(length(p), obj.NofPhases);
            for i = 1:obj.NofPhases
                rho(:, i) = State.Properties(['rho_', num2str(i)]).Value;
            end
            
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                for j=1:length(a)
                    Jnwp(a(j),a(j)) = Jnwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 2)*Inj(i).rho(j, 2);
                    Jwp(a(j),a(j)) = Jwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 1)*Inj(i).rho(j, 1);
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                for j=1:length(b)
                    Jnwp(b(j),b(j)) = Jnwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 2) .* rho(b(j), 2)...
                     - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * obj.drhodp(b(j),2) .* (Prod(i).p(j) - p(b(j)));                    
                    Jwp(b(j),b(j)) = Jwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 1) .* rho(b(j), 1)...
                     - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) .* obj.drhodp(b(j),1) .* (Prod(i).p(j) - p(b(j)));
                    
                    JwS(b(j),b(j)) = JwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* rho(b(j), 1) .* (Prod(i).p(j) - p(b(j))).*obj.dMob(b(j), 1);
                    JnwS(b(j),b(j)) = JnwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* rho(b(j), 2) .* (Prod(i).p(j) - p(b(j))).*obj.dMob(b(j), 2);
                end
            end
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, R, P)
            % Perform Average for ADM
            S_rest = R * Status.S;
            Status.S = P * S_rest;
            % Update other unknwons as well 
            %obj.UpdatePhaseCompositions(Status, FluidModel);
        end
        function CFL = ComputeCFLNumber(obj, ProductionSystem, DiscretizationModel, dt)
            N = DiscretizationModel.ReservoirGrid.N;      
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            P = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            
            % Copy values in local variables
            for i=1:obj.NofPhases
                P(:, i) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
                S(:, i) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(i)]).Value;
            end
            
            for i=1:obj.NofPhases
                [obj.Tph{i}, obj.Gph{i}] = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, obj.UpWind{i, 1}, obj.Mob(1:N,i), rho(:,i), obj.GravityModel.RhoInt{i, 1});
            end
             % Depths
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            % Source terms
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            ThroughPut = zeros(N, obj.NofPhases);
            Mass = zeros(N, obj.NofPhases);
            for i=1:obj.NofPhases
                % extract lower diagonals
                d = tril(obj.Tph{i}, -1);
                % extract upper diagonals
                u = triu(obj.Tph{i},  1);
                % extract main diagonal
                Diag = diag(obj.Tph{1});
                % Assemble
                D = diag(Diag + sum(u, 2)) + d; 
                U = diag(Diag + sum(d, 2)) + u;
                % extract lower diagonals
                gd = tril(obj.Gph{i}, -1);
                % extract upper diagonals
                gu = triu(obj.Gph{i},  1);
                % extract main diagonal
                DiagG = diag(obj.Gph{1});
                % Assemble
                GD = diag(DiagG + sum(gu, 2)) + gd; 
                GU = diag(DiagG + sum(gd, 2)) + gu;
                ThroughPut(:,i) = ...
                           - min(D * P(:,i), 0)...    % Convective term (take only incoming fluxes (negative))                
                           - min(U * P(:,i), 0)...    % Convective term (take only incoming fluxes (negative))     
                           - min(GD * depth, 0)...    % Gravity term (take only incoming fluxes (negative))                
                           - min(GU * depth, 0)...    % Gravity term (take only incoming fluxes (negative))
                           + max(q(:,i), 0);          % Wells (injectors)
                Mass(:,i) = rho(:,i) .* S(:,i) * pv;
            end
            Mass = max(Mass, 1e-10);
            ThroughPut(ThroughPut < 1e-10) = 0;
            Ratio = ThroughPut ./ Mass;
            CFL = dt * max(max(Ratio));
        end
        %% Methods for Sequential Coupling
        function ComputeTotalMobility(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.Mobt = sum(obj.Mob, 2);
        end
        function UpdateFractionalFlow(obj, ProductionSystem, FluidModel)
            obj.ComputeTotalMobility(ProductionSystem, FluidModel);
            rho = 0 * obj.Mob; 
            for i=1:obj.NofPhases
                rho(:,i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end
            obj.f = rho(:, 1) .* obj.Mob(:,1) ./ (rho(:, 1) .*obj.Mob(:,1) + rho(:, 2) .* obj.Mob(:, 2));
        end
        function dfdS(obj, ProductionSystem, FluidModel)
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            rho = 0 * obj.Mob; 
            for i=1:obj.NofPhases
                rho(:,i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end
            num = rho(:, 1) .* obj.Mob(:, 1);
            dnum = rho(:, 1) .* obj.dMob(:,1);
            den = sum(rho .* obj.Mob, 2);
            dden = sum(rho .* obj.dMob, 2);  
            obj.df = (dnum .* den - dden .* num) ./ den.^2;
            %obj.df = (obj.dMob(:,1) .* sum(obj.Mob, 2) - sum(obj.dMob, 2) .* obj.Mob(:,1)) ./ sum(obj.Mob, 2).^2;
        end
        function dfdsds = dfdSdS(obj, s, FluidModel)
           % f = Mob1 / (Mob1 + Mob2);
           % df = (dMob1 * (Mob1 + Mob2) - (dMob1 + dMob2)*Mob1 )/(Mob1 +
           % Mob2)^2
           % d(df)
           Mob = FluidModel.ComputePhaseMobilities(s);
           dMob = FluidModel.DMobDS(s);
           ddMob = FluidModel.DMobDSDS(s);
           A = dMob(:,1) .* sum(Mob, 2);
           dA = ddMob(:,1) .* sum(Mob, 2) +  sum(dMob, 2) .* dMob(:, 1);
           B = sum(dMob, 2) .* Mob(:,1);
           dB = sum(ddMob, 2) .* Mob(:,1) + sum(dMob, 2) .* ddMob(:,1);
           num = A - B;
           dnum = dA - dB;
           den = sum(Mob, 2).^2;
           dden = sum(Mob, 2) .* sum(dMob, 2);  
           
           dfdsds = (dnum .* den - dden .* num) ./ den.^2;
           
        end
        function Residual = BuildPressureResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            %%%% fix wells
            % Compute vector of qs
            qw = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            qf = zeros(DiscretizationModel.N, obj.NofPhases);
            if ProductionSystem.FracturesNetwork.Active
                qf = ComputeSourceTerms_frac_mat(obj, ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures, DiscretizationModel);
            end
            %% Reservoir residual
            Residual = zeros(DiscretizationModel.N, 1);
            % Transmissibility of reservoir
            for i=1:obj.NofPhases
                [obj.Tph{i, 1}, obj.Gph{i, 1}] = ...
                    obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, obj.UpWind{i, 1}, obj.Mob(1:DiscretizationModel.ReservoirGrid.N, i), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(i)]).Value, obj.GravityModel.RhoInt{i, 1});
            end
            Index.Start = 1;
            Index.End = DiscretizationModel.ReservoirGrid.N;
            phi = ProductionSystem.Reservoir.Por;
            Residual(1:DiscretizationModel.ReservoirGrid.N) =...
                obj.MediumPressureResidual(dt, phi, DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir.State, State0, qw, qf, Index, 0);
            if ProductionSystem.FracturesNetwork.Active
                %% Fractures Residuals
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = DiscretizationModel.Index_Local_to_Global(DiscretizationModel.ReservoirGrid.Nx, DiscretizationModel.ReservoirGrid.Ny, DiscretizationModel.ReservoirGrid.Nz, f, 1);
                    Index.End = DiscretizationModel.Index_Local_to_Global(DiscretizationModel.ReservoirGrid.Nx, DiscretizationModel.ReservoirGrid.Ny, DiscretizationModel.ReservoirGrid.Nz, f, DiscretizationModel.FracturesGrid.Grids(f).N);
                    % Transmissibility of fractures cells
                    for i=1:obj.NofPhases
                        [obj.Tph{i, 1+f}, obj.Gph{i, 1+f}] = ...
                            obj.TransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), obj.UpWind{i, 1+f}, obj.Mob(Index.Start:Index.End, i),...
                            ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(i)]).Value, obj.GravityModel.RhoInt{i, f+1});
                    end
                    Por = ProductionSystem.FracturesNetwork.Fractures(f).Por;
                    % Residual of Fracture Grid-cells
                    Residual(Index.Start:Index.End) = ...
                        obj.MediumPressureResidual(dt, phi, DiscretizationModel.FracturesGrid.Grids(f), ProductionSystem.FracturesNetwork.Fractures(f).State, State0, qw, qf, Index, f);
                end
            end
        end
        function Residual = MediumPressureResidual(obj, dt, phi, Grid, State, State0, qw, qf, Index, f)
            % Initialise local variables
            N = Grid.N;
            s_old = zeros(N, obj.NofPhases);
            rho_old = zeros(N, obj.NofPhases);
            P = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            pv = phi * Grid.Volume;
            
            % Copy values in local variables

            for i = 1:obj.NofPhases
                P(:, i) = State.Properties(['P_', num2str(i)]).Value;
                s_old(:, i) = State0.Properties(['S_', num2str(i)]).Value(Index.Start:Index.End);
                rho_old(:, i) = State0.Properties(['rho_', num2str(i)]).Value(Index.Start:Index.End);
                rho(:, i) = State.Properties(['rho_', num2str(i)]).Value;
                obj.Tph{i, f+1} = bsxfun(@rdivide, obj.Tph{i, f+1}, rho(:, i));
            end 
            depth = Grid.Depth;
            
            Residual = zeros(N, 1);
            for i=1:obj.NofPhases
            Residual(:) = ...
                    Residual(:) - ...
                    pv/dt * (rho_old(:, i) .* s_old(:, i) ./ rho(:, i))...
                    + obj.Tph{i, f+1} * P(:, i)...
                    - qw(:, i) ./ rho(:, i)...
                    - qf(:, i) ./ rho(:, i);
            end
            Residual = Residual + pv/dt;
%             
%             Residual(:) = Residual(:) + pv/dt; 
%             % without dividing 'Mattteo'
%             Residual = zeros(N, 1);
%             for i=1:obj.NofPhases
%             Residual(:) = ...
%                    Residual(:) + ...
%                    pv/dt * (rho(:, i).*s_old(:, i) - rho_old(:, i) .* s_old(:, i))...
%                    + obj.Tph{i} * P(:, i)...
%                    - qw(:, i)...
%                    - qf(:, i);
%             end
        end
        function A = BuildPressureMatrix(obj, ProductionSystem, DiscretizationModel, dt, State0)
            %% 1. Reservoir            
            A = obj.Build_Medium_PressureMatrix(ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, 1, DiscretizationModel.ReservoirGrid.N, 0);
            %% 3. Add wells
            rho = zeros(DiscretizationModel.ReservoirGrid.N, obj.NofPhases);
            for i=1:obj.NofPhases
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end
            W = obj.AddWellsToPressureSystem(DiscretizationModel.ReservoirGrid.N, ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K(:,1), rho);
            A = A + W;
           
            %% 2. Fractures
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, DiscretizationModel.FracturesGrid.Grids(f).N);
                Af = obj.Build_Medium_PressureMatrix(ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Start, End, f);
                A = blkdiag(A, Af);
            end
            %% 3.  Add matrix-frac and frac-frac connections
            for If1_Local = 1 : length(DiscretizationModel.CrossConnections)
                for i = 1:obj.NofPhases
                    % |x    |  x  |
                    % |  M  |  MF |
                    % |  x  |  x  |
                    % |-----|-----|----
                    % | x  x|s(x) |
                    % |  FM |  F  |
                    % |     |     |
                    
                    % matrix-frac
                    % Fill in FM entries
                    indices_m = DiscretizationModel.CrossConnections(If1_Local).Cells( DiscretizationModel.CrossConnections(If1_Local).Cells <= DiscretizationModel.ReservoirGrid.N );
                    If1_Global = DiscretizationModel.ReservoirGrid.N+If1_Local; % Global index of this fracture cells;
                    Index_frac1_Local = DiscretizationModel.Index_Global_to_Local(If1_Global);
                    % Fill in FM entries
                    A(If1_Global, indices_m) = A(If1_Global, indices_m) - ...
                        DiscretizationModel.CrossConnections(If1_Local).T_Geo(1:length(indices_m))'...
                        .* (obj.Mob(indices_m,i)' .* ProductionSystem.Reservoir.State.Properties(['rho_',num2str(i)]).Value(indices_m)' );
                    % .*(obj.Mob(If1_Global,i) .* ProductionSystem.FracturesNetwork.Fractures(Index_frac1_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac1_Local.g) ...
                    % Fill in F diagonl entries
                    A(If1_Global, If1_Global) = A(If1_Global, If1_Global) + ...
                        sum( DiscretizationModel.CrossConnections(If1_Local).T_Geo(1:length(indices_m)) .* obj.Mob(If1_Global,i) .* ProductionSystem.Reservoir.State.Properties(['rho_',num2str(i)]).Value(indices_m) );
                    % Fill in MF entries
                    A(indices_m, If1_Global) = A(If1_Global) - ...
                        DiscretizationModel.CrossConnections(If1_Local).T_Geo(1:length(indices_m)) .* obj.Mob(If1_Global,i) .* ProductionSystem.Reservoir.State.Properties(['rho_',num2str(i)]).Value(indices_m);
                    % Fill in M diagonal entries
                    A(sub2ind(size(A), indices_m, indices_m)) = A(sub2ind(size(A), indices_m, indices_m)) + ...
                        DiscretizationModel.CrossConnections(If1_Local).T_Geo(1:length(indices_m)) .* obj.Mob(If1_Global,i) .* ProductionSystem.Reservoir.State.Properties(['rho_',num2str(i)]).Value(indices_m);
                    
                    % frac-frac
                    indices_f = DiscretizationModel.CrossConnections(If1_Local).Cells( DiscretizationModel.CrossConnections(If1_Local).Cells > DiscretizationModel.ReservoirGrid.N );
                    if ~isempty(indices_f)
                        for n = 1:length(indices_f)
                            If2_Global = indices_f(n); % Global indices of the other fractures' cells if any
                            If2_Local = If2_Global - DiscretizationModel.ReservoirGrid.N;
                            T_Geo_Half1 = DiscretizationModel.CrossConnections(If1_Local      ).T_Geo( length(indices_m)+n );
                            T_Geo_Half2 = DiscretizationModel.CrossConnections(If2_Local).T_Geo( DiscretizationModel.CrossConnections(If2_Local).Cells==If1_Global );
                            T_Geo_Harmo = (T_Geo_Half1 * T_Geo_Half2) / (T_Geo_Half1 + T_Geo_Half2);
                            Index_frac2_Local = DiscretizationModel.Index_Global_to_Local(If2_Global);
                            % Fill in FG entries
                            A(If1_Global, If2_Global) = A(If1_Global, If2_Global) - ...
                                T_Geo_Harmo * obj.Mob(If1_Global,i) * ProductionSystem.FracturesNetwork.Fractures(Index_frac1_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac1_Local.g);
                            % Fill in F diagonal entries
                            A(If1_Global, If1_Global) = A(If1_Global, If1_Global) + ...
                                T_Geo_Harmo * obj.Mob(If1_Global,i) * ProductionSystem.FracturesNetwork.Fractures(Index_frac1_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac1_Local.g);
                            % Fill in GF diagonal entries
                            A(If2_Global, If1_Global) = A(If2_Global, If1_Global) - ...
                                T_Geo_Harmo * obj.Mob(If1_Global,i) * ProductionSystem.FracturesNetwork.Fractures(Index_frac1_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac1_Local.g);
                            % Fill in G diagonal entries
                            A(If2_Global, If2_Global) = A(If2_Global, If2_Global) + ...
                                T_Geo_Harmo * obj.Mob(If1_Global,i) * ProductionSystem.FracturesNetwork.Fractures(Index_frac1_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac1_Local.g);
                        end
                    end
                end
            end
        end
        function J = Build_Medium_PressureMatrix(obj, Medium, Grid, dt, State0, Start, End, f)
            %% 0. Initialise local variables
            N = Grid.N;
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            
            pv = Medium.Por*Grid.Volume;
            s_old = zeros(N, obj.NofPhases);
            rho_old = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            
            % Copy values in local variables
            for i=1:obj.NofPhases
                s_old(:, i) = State0.Properties(['S_', num2str(i)]).Value;
                rho_old(:, i) = State0.Properties(['rho_', num2str(i)]).Value;
                rho(:, i) = Medium.State.Properties(['rho_', num2str(i)]).Value;
            end 
            
            A = sparse(N,N);
            C = sparse(N,N);
            for i=1:obj.NofPhases
                %% 1. Accummulation term
                vec = - pv/dt * (- rho_old(Start:End, i) .* s_old(Start:End, i) .* obj.drhodp(Start:End,i) ./ rho(:, i).^2) ;
                % Alternative
                %vec = pv/dt * obj.drhodp(Start:End,i);
                C = C + spdiags(vec, 0, N, N);
                
                %% 2. Convective term 
                % a. transmissibility matrix
                A = A + obj.Tph{i, f+1};
                % b. Derivative of transmissibility with respect to p
                dMupx = obj.UpWind{i, f+1}.x*(obj.Mob(Start:End, i) .* obj.drhodp(Start:End,i));
                dMupy = obj.UpWind{i, f+1}.y*(obj.Mob(Start:End, i) .* obj.drhodp(Start:End,i));
                dMupz = obj.UpWind{i, f+1}.z*(obj.Mob(Start:End, i) .* obj.drhodp(Start:End,i));
                
                vecX1 = min(reshape(obj.U{i, f+1}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{i, f+1}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{i, f+1}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{i, f+1}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{i, f+1}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{i, f+1}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz; 
                maindiag = vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1;
                
                DiagVecs = [-vecZ2, -vecY2, -vecX2, maindiag, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                dTdp = spdiags(DiagVecs, DiagIndx, N, N);
%                 
                % multiply row i by rho(i)
                dTdp = bsxfun(@times, dTdp, rho(:, i));
                
%                 dMupx = obj.UpWind{i, f+1}.x*(obj.Mob(:, i) .* rho(:,i));
%                 dMupy = obj.UpWind{i, f+1}.y*(obj.Mob(:, i) .* rho(:,i));
%                 dMupz = obj.UpWind{i, f+1}.z*(obj.Mob(:, i) .* rho(:,i));
%                 
%                 vecX1 = min(reshape(obj.U{i, f+1}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
%                 vecX2 = max(reshape(obj.U{i, f+1}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
%                 vecY1 = min(reshape(obj.U{i, f+1}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
%                 vecY2 = max(reshape(obj.U{i, f+1}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
%                 vecZ1 = min(reshape(obj.U{i, f+1}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
%                 vecZ2 = max(reshape(obj.U{i, f+1}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz; 
%                 maindiag = vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1;
%                 vec = - maindiag .* obj.drhodp(:,i); 
%                 dTdp = dTdp + spdiags(vec, 0, N, N);
                 
                dTdp = bsxfun(@rdivide, dTdp, rho(:, i).^2);
%                  
                 A = A + dTdp; 
            end            
            %% Jacobian: Put them together
            J = A + C;
            
        end       
        function W = AddWellsToPressureSystem(obj, N, State, Wells, K, rho)
            %% Add Wells in residual form
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            dq = zeros(N, obj.NofPhases);
            %Injectors
            for i=1:length(Inj)
                c = Inj(i).Cells;
                [dq(c, :), ~] = Inj(i).dQPhasesdPdS(State, K, obj.NofPhases);
            end
            
            %Producers
            for i=1:length(Prod)
                c = Prod(i).Cells;
                [dq(c, :), ~] = Prod(i).dQPhasesdPdS(State, K, obj.Mob, obj.dMob, obj.drhodp, obj.NofPhases);
            end
            q = sparse(obj.ComputeSourceTerms(N, Wells));
            dqt = sum((dq .* rho - obj.drhodp .* q)./rho.^2, 2);
            % Alternative
            %dqt = sum(dq, 2);
            
            W = -spdiags(dqt, 0, N, N);
            
        end
        function UpdatePressure(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            %% 1. Update matrix pressure and densities
            % Update Pressure
            Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(obj.NofPhases)]);
            Pm.update(delta(1:DiscretizationModel.ReservoirGrid.N));
            % Update Phase Densities
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            % Update total density
            FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State);
            if obj.NofPhases > 1
            % Update Pc
                FluidModel.ComputePc(ProductionSystem.Reservoir.State);
            end
            
            %% 2. Update fractures pressure and densities
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    % Update Pressure
                    index1 = DiscretizationModel.ReservoirGrid.N + sum(DiscretizationModel.FracturesGrid.N(1:f-1)) + 1;
                    index2 = index1 - 1 + DiscretizationModel.FracturesGrid.N(f);
                    Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(obj.NofPhases)]);
                    Pf.update(delta(index1:index2));
                    % Update Phase Densities
                    FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    % Update total density
                    FluidModel.ComputeTotalDensity(ProductionSystem.FracturesNetwork.Fractures(f).State);
                end
            end
        end
        function ComputeTotalFluxes(obj, ProductionSystem, DiscretizationModel)
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            
            obj.Utot.x = zeros(Nx+1, Ny, Nz); 
            obj.Utot.y = zeros(Nx, Ny+1, Nz); 
            obj.Utot.z = zeros(Nx, Ny, Nz+1);
            
            N = Nx*Ny*Nz;
            rho = zeros(N, obj.NofPhases);
            for i=1:obj.NofPhases
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end 
            for ph = 1:obj.NofPhases
                obj.Utot.x(2:Nx+1,:,:) = obj.Utot.x(2:Nx+1,:,:) + obj.U{ph, 1}.x(2:Nx+1,:,:) .* reshape(obj.UpWind{ph,1}.x *  (rho(:, ph) .* obj.Mob(:, ph)), Nx, Ny, Nz); %- Ucap.x(2:Nx,:);
                obj.Utot.y(:,2:Ny+1,:) = obj.Utot.y(:,2:Ny+1,:) + obj.U{ph, 1}.y(:,2:Ny+1,:) .* reshape(obj.UpWind{ph,1}.y *  (rho(:, ph) .* obj.Mob(:, ph)), Nx, Ny, Nz); %- Ucap.y(:,2:Ny);
                obj.Utot.z(:,:,2:Nz+1) = obj.Utot.z(:,:,2:Nz+1) + obj.U{ph, 1}.z(:,:,2:Nz+1) .* reshape(obj.UpWind{ph,1}.z *  (rho(:, ph) .* obj.Mob(:, ph)), Nx, Ny, Nz);  %- Ucap.y(:,2:Ny);
            end
            
            % Wells total fluxes
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            obj.Qwells = sum(q, 2);
        end
        function conservative = CheckMassConservation(obj, Grid)
            %Checks mass balance in all cells
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            conservative = 1;
            maxUx = max(max(max(obj.Utot.x)));
            maxUy = max(max(max(obj.Utot.y)));
            maxUz = max(max(max(obj.Utot.z)));
            maxU = max([maxUx, maxUy, maxUz]);
            qWells = reshape(obj.Qwells, Nx, Ny, Nz);
            for k=1:Nz
                for j=1:Ny
                    for i=1:Nx
                        Accum = obj.Utot.x(i,j,k) - obj.Utot.x(i+1,j,k) + obj.Utot.y(i,j,k) - obj.Utot.y(i,j+1,k) + obj.Utot.z(i,j,k) - obj.Utot.z(i,j,k+1) + qWells(i,j,k);
                        if (abs(Accum/maxU) > 10^(-5))
                            conservative = 0;
                        end
                    end
                end
            end
        end
        function ViscousMatrix(obj, Grid)
            %Builds Upwind Flux matrix
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;                                   
            q = min(obj.Qwells, 0);                        
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(obj.Utot.x, 0); 
            Yneg = min(obj.Utot.y, 0);
            Zneg = min(obj.Utot.z, 0);
            % make them vectors 
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(obj.Utot.x, 0); 
            Ypos = max(obj.Utot.y, 0); 
            Zpos = max(obj.Utot.z, 0);
            % make them vectors
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
            % Assemble matrix
            DiagVecs = [z2, y2, x2, q+x1-x2+y1-y2+z1-z2, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny]; % diagonal index
            obj.V = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function Residual = BuildTransportResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Initialise local objects
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            s = ProductionSystem.Reservoir.State.Properties('S_1').Value;
            s_old = State0.Properties('S_1').Value;      
            rho = ProductionSystem.Reservoir.State.Properties('rho_1').Value;
            
            % viscous fluxes matrix
            obj.ViscousMatrix(DiscretizationModel.ReservoirGrid);      
            
            % Compute residual
            Residual = pv/dt .* rho .* (s - s_old)  - max(obj.Qwells, 0) - obj.V * obj.f;
        end
        function Jacobian = BuildTransportJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            % Build Transport Jacobian
            N = DiscretizationModel.ReservoirGrid.N;
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            rho = ProductionSystem.Reservoir.State.Properties('rho_1').Value;
            
            D = spdiags(pv/dt*rho, 0, N, N);
            Jacobian = D - obj.V * spdiags(obj.df,0,N,N); %+ CapJac;
        end
        function delta = UpdateSaturation(obj, ProductionSystem, delta, FluidModel, DiscretizationModel)
            s_old = ProductionSystem.Reservoir.State.Properties('S_1').Value; 
            
            % Update
            snew = s_old + delta;
            
            % Remove values that are not physical
            snew = max(snew, 0);
            snew = min(snew, 1);
            
            % FLUX CORRECTION - PATRICK
            Ddf_old = obj.dfdSdS(s_old, FluidModel);
            Ddf = obj.dfdSdS(snew, FluidModel);           
            snew = snew.*(Ddf.*Ddf_old >= 0) + 0.5*(snew + s_old).*(Ddf.*Ddf_old<0);
            delta = snew-s_old;
            
            % This is the actual update
            Nm = DiscretizationModel.ReservoirGrid.N;
            DeltaLast = zeros(Nm, 1);
            for ph = 1:obj.NofPhases-1
                Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(ph)]);
                DeltaS = delta(1:Nm);
                Sm.update(DeltaS);
                % Remove values that are not physical
                Sm.Value = max(Sm.Value, 0);
                Sm.Value = min(Sm.Value, 1);
                DeltaLast = DeltaLast + DeltaS;
            end
            Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(obj.NofPhases)]);
            Sm.update(-DeltaLast);
        end
        function delta = UpdateSaturationExplicitly(obj, ProductionSystem, DiscretizationModel, dt)
            % 0. Initialise
            pv = ProductionSystem.Reservoir.Por .* DiscretizationModel.ReservoirGrid.Volume;
            N = DiscretizationModel.ReservoirGrid.N;
            
            % 1. Solve
            T = spdiags(dt/pv*ones(N,1),0,N,N);    % dt/pv * Cell Fluxes and producer
            B = T * obj.V;
            injector = max(obj.Qwells,0) .* dt/pv;  % injection flux * dt/pv
            
            S = ProductionSystem.Reservoir.State.Properties('S_1');
            s_old = S.Value;
            
            delta = S.Value - s_old + (B * obj.f + injector);
            
            % Now update values
            Nm = DiscretizationModel.ReservoirGrid.N;
            DeltaLast = zeros(Nm, 1);
            for ph = 1:obj.NofPhases-1
                Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(ph)]);
                DeltaS = delta(1:Nm);
                Sm.update(DeltaS);
                % Remove values that are not physical
                Sm.Value = max(Sm.Value, 0);
                Sm.Value = min(Sm.Value, 1);
                DeltaLast = DeltaLast + DeltaS;
            end
            Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(obj.NofPhases)]);
            Sm.update(-DeltaLast);
            % Remove values that are not physical
            Sm.Value = max(Sm.Value, 0);
            Sm.Value = min(Sm.Value, 1);
        end
    end
end