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
        Mobt
    end
    methods
        function obj = Immiscible_formulation()
            obj@formulation();
            obj.Tph = cell(2,1);
            obj.Gph = cell(2,1);
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
            Jacobian = BuildMediumJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0);
            
            % Jacobian of the fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Nf = DiscretizationModel.FracturesGrid.N;
                Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                Index.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, DiscretizationModel.FracturesGrid.Grids(f).N);
                Jacobian_frac_f = BuildMediumJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f);
                Jacobian  = blkdiag( Jacobian , Jacobian_frac_f );
            end
            
            % Effects of frac-matrix & frac-frac connections in Jacobian
            for If1_Local = 1 : length(DiscretizationModel.CrossConnections)
                DiscretizationModel.CrossConnections(If1_Local).CrossConnectionsUpWind(ProductionSystem, DiscretizationModel, If1_Local, obj);
                Cells = DiscretizationModel.CrossConnections(If1_Local).Cells;
                T_Geo = DiscretizationModel.CrossConnections(If1_Local).T_Geo;
                UpWind = DiscretizationModel.CrossConnections(If1_Local).UpWind;
                for i = 1:obj.NofPhases
                    % Loop over equations
                    If1_Global = If1_Local + Nm;
                    Index_frac1_Local = DiscretizationModel.Index_Global_to_Local(If1_Global);
                    rho_f1 = Fractures(Index_frac1_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac1_Local.g);
                    Pf1 = Fractures(Index_frac1_Local.f).State.Properties(['P_',num2str(i)]).Value(Index_frac1_Local.g);

                    % frac-matrix
                    % Pressure Blocks
                    indices_m = Cells( Cells <= Nm );
                    rho_m = Reservoir.State.Properties(['rho_',num2str(i)]).Value(indices_m);
                    Pm = Reservoir.State.Properties(['P_',num2str(i)]).Value(indices_m);
                    Jp_frac_mat = - T_Geo(1:length(indices_m)) .* (...
                         UpWind(1:length(indices_m),i).*obj.Mob(indices_m,i) .*( rho_m  + obj.drhodp(indices_m,i) .*(Pm-Pf1) ) + ...
                        ~UpWind(1:length(indices_m),i).*obj.Mob(If1_Global,i).*( rho_f1 + obj.drhodp(If1_Global,i).*(Pm-Pf1) ) );
                    
                    % Adding the Jp_frac_mat of phase i to Jacobian
                    % adding to Jp_FM
                    row_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                    col_ind_Jp = indices_m;
                    Jacobian(row_ind_Jp, col_ind_Jp) = Jp_frac_mat';
                    % adding to Jp_FF
                    Jacobian(row_ind_Jp, row_ind_Jp - (i-1)*Nf(Index_frac1_Local.f)) = ...
                    Jacobian(row_ind_Jp, row_ind_Jp - (i-1)*Nf(Index_frac1_Local.f)) - sum(Jp_frac_mat);
                    % adding to Jp_MF
                    col_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + Index_frac1_Local.g;
                    row_ind_Jp = (i-1)*Nm + indices_m;
                    Jacobian(row_ind_Jp, col_ind_Jp) = Jp_frac_mat;
                    % adding to Jp_MM
                    Jacobian(sub2ind([Nt*obj.NofPhases-(i-1)*Nm , Nm], indices_m, indices_m) + (i-1).*Nm.*indices_m) = ...
                    Jacobian(sub2ind([Nt*obj.NofPhases-(i-1)*Nm , Nm], indices_m, indices_m) + (i-1).*Nm.*indices_m) - Jp_frac_mat;
  
                    % Saturation Blocks
                    Js_mat_frac = T_Geo(1:length(indices_m)) .*  UpWind(1:length(indices_m),i).*(Pf1-Pm).*( rho_m .* obj.dMob(indices_m ,i) );
                    Js_frac_mat = T_Geo(1:length(indices_m)) .* ~UpWind(1:length(indices_m),i).*(Pf1-Pm).*( rho_f1.* obj.dMob(If1_Global,i) );
                    % Adding the Js_frac_mat of phase i to matrix columns
                    % adding to Js_FM
                    row_ind_Js = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                    col_ind_Js = indices_m + Nm;
                    Jacobian(row_ind_Js, col_ind_Js) = Js_mat_frac;
                    % adding to Js_MM
                    Jacobian(sub2ind([Nt*obj.NofPhases-(i-1)*Nm , Nm], indices_m, indices_m) + obj.NofPhases*Nm*Nt + (i-1).*Nm.*indices_m) = ...
                    Jacobian(sub2ind([Nt*obj.NofPhases-(i-1)*Nm , Nm], indices_m, indices_m) + obj.NofPhases*Nm*Nt + (i-1).*Nm.*indices_m) - Js_mat_frac;
                    % Adding the Js_frac_mat of phase i to fractures columns
                    % adding to Js_MF
                    row_ind_Js = (i-1)*Nm + indices_m;
                    col_ind_Js = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + Nf(Index_frac1_Local.f) + Index_frac1_Local.g;
                    Jacobian(row_ind_Js, col_ind_Js) = -Js_frac_mat;
                    % adding to Js_FF
                    row_ind_Js = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f)) + Index_frac1_Local.g;
                    Jacobian(row_ind_Js, row_ind_Js + (obj.NofPhases-i)*Nf(Index_frac1_Local.f)) = ...
                    Jacobian(row_ind_Js, row_ind_Js + (obj.NofPhases-i)*Nf(Index_frac1_Local.f)) + sum(Js_frac_mat);
                    
                    % frac-frac
                    indices_f = DiscretizationModel.CrossConnections(If1_Local).Cells( DiscretizationModel.CrossConnections(If1_Local).Cells > Nm );
                    if ~isempty(indices_f)
                        for n = 1:length(indices_f)
                            If2_Global = indices_f(n); % Global indices of the other fractures' cells if any
                            If2_Local = If2_Global - Nm;
                            Index_frac2_Local = DiscretizationModel.Index_Global_to_Local(If2_Global);
                            T_Geo_Half1 = DiscretizationModel.CrossConnections(If1_Local).T_Geo( length(indices_m)+n );
                            rho_f2 = Fractures(Index_frac2_Local.f).State.Properties(['rho_',num2str(i)]).Value(Index_frac2_Local.g);
                            Pf2 = Fractures(Index_frac2_Local.f).State.Properties(['P_',num2str(i)]).Value(Index_frac2_Local.g);
                            T_Geo_Half2 = DiscretizationModel.CrossConnections(If2_Local).T_Geo( DiscretizationModel.CrossConnections(If2_Local).Cells==If1_Global );
                            T_Geo_Harmo = (T_Geo_Half1 * T_Geo_Half2) / (T_Geo_Half1 + T_Geo_Half2);
                            
                            % Pressure Blocks
                            Jp_frac_frac = - T_Geo_Harmo * (...
                                UpWind(length(indices_m)+n,i) * obj.Mob(If1_Global,i) * ( rho_f2 + obj.drhodp(If2_Global,i)*(Pf2-Pf1) ) + ...
                               ~UpWind(length(indices_m)+n,i) * obj.Mob(If2_Global,i) * ( rho_f1 + obj.drhodp(If1_Global,i)*(Pf2-Pf1) ) );
                            % Adding the Jp_frac_frac of phase i to Jacobian
                            % adding to Jp_F1F2
                            row_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f-1)) + Index_frac1_Local.g;
                            col_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) +                                   + Index_frac2_Local.g;
                            Jacobian(row_ind_Jp, col_ind_Jp) = Jp_frac_frac;
                            % adding to Jp_F1F1
                            row_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f-1)) + Index_frac1_Local.g;
                            col_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) +                                   + Index_frac1_Local.g;
                            Jacobian(row_ind_Jp, col_ind_Jp) = Jacobian(row_ind_Jp, col_ind_Jp) - Jp_frac_frac;
                            % adding to Jp_F2F1
                            row_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) + (i-1)*(Nf(Index_frac2_Local.f-1)) + Index_frac2_Local.g;
                            col_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) +                                   + Index_frac1_Local.g;
                            Jacobian(row_ind_Jp, col_ind_Jp) = Jp_frac_frac;
                            % adding to Jp_F2F2
                            row_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) + (i-1)*(Nf(Index_frac2_Local.f-1)) + Index_frac2_Local.g;
                            col_ind_Jp = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) +                                   + Index_frac2_Local.g;
                            Jacobian(row_ind_Jp, col_ind_Jp) = Jacobian(row_ind_Jp, col_ind_Jp) - Jp_frac_frac;
                           
                            % Saturation Blocks
                            Js_frac2_frac1 = T_Geo_Harmo* UpWind(length(indices_m)+n,i)*(Pf1-Pf2)*(rho_f2*obj.dMob(If2_Global,i));
                            Js_frac1_frac2 = T_Geo_Harmo*~UpWind(length(indices_m)+n,i)*(Pf1-Pf2)*(rho_f1*obj.dMob(If1_Global,i));
                            % Adding the Js_frac2_frac1 of phase i to frac2 columns
                            % adding to Js_F1F2
                            row_ind_Js = obj.NofPhases*( Nm + sum(Nf(1:Index_frac1_Local.f-1)) ) + (i-1)*(Nf(Index_frac1_Local.f-1)) + Index_frac1_Local.g;
                            col_ind_Js = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) +    1 *(Nf(Index_frac2_Local.f  )) + Index_frac2_Local.g;
                            Jacobian(row_ind_Js, col_ind_Js) = Js_frac2_frac1;
                            % adding to Js_F2F2
                            row_ind_Js = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) + (i-1)*(Nf(Index_frac2_Local.f-1)) + Index_frac2_Local.g;
                            col_ind_Js = obj.NofPhases*( Nm + sum(Nf(1:Index_frac2_Local.f-1)) ) +     1*(Nf(Index_frac2_Local.f  )) + Index_frac2_Local.g;
                            Jacobian(row_ind_Js, col_ind_Js) = Jacobian(row_ind_Js, col_ind_Js) - Js_frac2_frac1;
                            % Adding the Js_frac1_frac2 of phase i to frac2 columns
                            % adding to Js_F2F1
                            
                        end
                    end
                end
            end
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
                    Sm.update(delta( ph*Nm+1 : (ph+1)*Nm ));
                    % Remove values that are not physical
                    Sm.Value = max(Sm.Value, 0);
                    Sm.Value = min(Sm.Value, 1);
                    DeltaLast = DeltaLast + delta( ph*Nm+1 : (ph+1)*Nm );
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
                        Pf.update(delta( obj.NofPhases*Nm+sum(Nf(1:f-1))+1 : obj.NofPhases*Nm+sum(Nf(1:f)) ));
                        DeltaLast = zeros(Nf(f), 1);
                        for ph = 1:obj.NofPhases-1
                            Sf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['S_', num2str(ph)]);
                            Sf.update(delta( obj.NofPhases*Nm+sum(Nf(1:f-1))+ Nf(f) + 1 : obj.NofPhases*Nm+sum(Nf(1:f)) + Nf(f)));
                            Sf.Value = max(Sf.Value, 0);
                            Sf.Value = min(Sf.Value, 1);
                            DeltaLast = DeltaLast + delta( obj.NofPhases*Nm+sum(Nf(1:f-1))+ Nf(f)+ 1 : obj.NofPhases*Nm+sum(Nf(1:f)) + Nf(f));
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
            CFL = 0;
        end
        %% Methods for Sequential Coupling
        function ComputeTotalMobility(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.Mobt = sum(obj.Mob, 2);
        end
        function UpdateFractionalFlow(obj, ProductionSystem, FluidModel)
            obj.ComputeTotalMobility(ProductionSystem, FluidModel);
            obj.f = obj.Mob(:,1) ./ obj.Mobt;
        end
        function dfdS(obj, ProductionSystem, FluidModel)
            dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.S);
            obj.df = (dMob(:,1) .* sum(obj.Mob, 2) - sum(dMob, 2) .* obj.Mob(:,1)) ./ sum(obj.Mob, 2).^2;
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
            pvdt = ProductionSystem.Reservoir.Por .* DiscretizationModel.ReservoirGrid.Volume ./ dt;
            Residual(1:DiscretizationModel.ReservoirGrid.N) =...
                obj.MediumPressureResidual(pvdt, DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir.State, State0, qw, qf, Index, 0);
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
                    pvdt = ProductionSystem.FracturesNetwork.Fractures(f).Por .* DiscretizationModel.FracturesGrid.Grids(f).Volume ./ dt;
                    % Residual of Fracture Grid-cells
                    Residual(Index.Start:Index.End) = ...
                        obj.MediumPressureResidual(pvdt, DiscretizationModel.FracturesGrid.Grids(f), ProductionSystem.FracturesNetwork.Fractures(f).State, State0, qw, qf, Index, f);
                end
            end
        end
        function Residual = MediumPressureResidual(obj, pvdt, Grid, State, State0, qw, qf, Index, f)
            % Initialise local variables
            N = Grid.N;
            s_old = zeros(N, obj.NofPhases);
            rho_old = zeros(N, obj.NofPhases);
            P = zeros(N, obj.NofPhases);
            s = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            
            % Copy values in local variables
            for i = 1:obj.NofPhases
                P(:, i) = State.Properties(['P_', num2str(i)]).Value;
                s_old(:, i) = State0.Properties(['S_', num2str(i)]).Value(Index.Start:Index.End);
                rho_old(:, i) = State0.Properties(['rho_', num2str(i)]).Value(Index.Start:Index.End);
                s(:, i) = State.Properties(['S_', num2str(i)]).Value;
                rho(:, i) = State.Properties(['rho_', num2str(i)]).Value;
            end 
            depth = Grid.Depth;
            
            % Accumulation Term
            AS = speye(N)*pvdt;
            
            %% RESIDUAL
            Residual = zeros(N,1);
            for i=1:obj.NofPhases
                Residual(:)  = Residual (:) + AS*(rho(:,i) .* s(:,i) - rho_old(:,i) .* s_old(:,i))...
                               + obj.Tph{i, f+1} * P(:, i)...
                               - qw(Index.Start:Index.End,i)...
                               - qf(Index.Start:Index.End,i);
            end
        end
        function A = BuildPressureMatrix(obj, ProductionSystem, DiscretizationModel, dt)
            %% 1. Reservoir
            A = obj.Tph{1, 1};
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            for i=2:obj.NofPhases
                A = A + obj.Tph{i, 1}; 
            end
            
            % Compressibility in the accumulation term (only for matrix)
            Index.Start = 1;
            Index.End = DiscretizationModel.ReservoirGrid.N;
            vec = pv .* obj.drhodp(Index.Start:Index.End) / dt;
            A = A + diag(vec);
            
            A = obj.AddWellsToPressureSystem(A, ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K(:,1));
            %% 2. Fractures
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                Af = obj.Tph{1 , f+1};
                for i=2:obj.NofPhases
                    Af = Af + obj.Tph{i, f+1};
                end
                pv = ProductionSystem.FracturesNetwork.Fractures(f).Por*DiscretizationModel.FracturesGrid.Grids(f).Volume; 
                Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, DiscretizationModel.FracturesGrid.Grids(f).N);
                vec = pv .* obj.drhodp(Start:End) / dt;
                Af = Af + diag(vec);
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
        function A = AddWellsToPressureSystem(obj, A, State, Wells, K)
            %% Add Wells in residual form
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                for j=1:length(a)
                    for phase=1:obj.NofPhases
                        A(a(j),a(j)) = A(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, phase)*Inj(i).rho(j, phase);
                    end
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                for j=1:length(b)
                    for phase=obj.NofPhases
                        A(b(j),b(j)) = A(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), phase) .* State.Properties(['rho_', num2str(phase)]).Value(b(j))...
                                       - Prod(i).PI * K(b(j)) * obj.Mob(b(j), phase) * obj.drhodp(b(j), phase) .* (Prod(i).p(j) - State.Properties(['P_', num2str(obj.NofPhases)]).Value(b(j)));                    
                    end
                end
            end
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
        function [Utot, Qwells] = ComputeTotalFluxes(obj, ProductionSystem, DiscretizationModel)
            Utot.x(2:Nx+1,:,:) = obj.U(2).x(2:Nx+1,:,:) .* reshape(obj.UpWind.x *  obj.Mobt, Nx, Ny, Nz); %- Ucap.x(2:Nx,:);
            Utot.y(:,2:Ny+1,:) = obj.U(2).y(:,2:Ny+1,:) .* reshape(obj.UpWind.y *  obj.Mobt, Nx, Ny, Nz); %- Ucap.y(:,2:Ny);
            Utot.z(:,:,2:Nz+1) = obj.U(2).z(:,:,2:Nz+1) .* reshape(obj.UpWind.z *  obj.Mobt, Nx, Ny, Nz);  %- Ucap.y(:,2:Ny);
            
            % Wells total fluxes
            Qwells = ProductionSystem.Wells.TotalFluxes(ProductionSystem.Reservoir, obj.Mobt);
        end
        function conservative = CheckMassConservation(obj, Grid, Utot, Qwells)
            %Checks mass balance in all cells
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            conservative = 1;
            maxUx = max(max(max(obj.Utot.x)));
            maxUy = max(max(max(obj.Utot.y)));
            maxUz = max(max(max(obj.Utot.z)));
            maxU = max([maxUx, maxUy, maxUz]);
            qWells = reshape(Qwells, Nx, Ny, Nz);
            for k=1:Nz
                for j=1:Ny
                    for i=1:Nx
                        Accum = Utot.x(i,j,k) - Utot.x(i+1,j,k) + Utot.y(i,j,k) - Utot.y(i,j+1,k) + Utot.z(i,j,k) - Utot.z(i,j,k+1) + qWells(i,j,k);
                        if (abs(Accum/maxU) > 10^(-5))
                            conservative = 0;
                        end
                    end
                end
            end
        end
        function ViscousMatrix(obj, Grid, Utot)
            %Builds Upwind Flux matrix
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;                                   
            q = min(obj.Qwells, 0);                        
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(Utot.x, 0); 
            Yneg = min(Utot.y, 0);
            Zneg = min(Utot.z, 0);
            % make them vectors 
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(Utot.x, 0); 
            Ypos = max(Utot.y, 0); 
            Zpos = max(Utot.z, 0);
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
            s = ProductionSystem.Reservoir.State.S;
            s_old = State0.S;      
            
            % Compute residual
            Residual = pv/dt * (s - s_old)  - max(obj.Qwells, 0) - obj.V * obj.f;
        end
        function Jacobian = BuildTransportJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            % Build Transport Jacobian
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            N = DiscretizationModel.ReservoirGrid.N;
            D = spdiags(pv/dt*ones(N,1),0,N,N);
            Jacobian = D - obj.V * spdiags(obj.df,0,N,N); %+ CapJac;
        end
        function UpdateSaturation(obj, State, delta, FluidModel)
            State.S = State.S + delta;
            State.S = min(State.S, 1);
            State.S = max(State.S, FluidModel.Phases(1).sr);
        end
        function UpdateSaturationExplicitly(obj, ProductionSystem, DiscretizationModel, dt)
            % 0. Initialise
            pv = ProductionSystem.Reservoir.Por .* DiscretizationModel.ReservoirGrid.Volume;
            N = DiscretizationModel.ReservoirGrid.N;
            
            % 1. Solve
            T = spdiags(dt/pv*ones(N,1),0,N,N);    % dt/pv * Cell Fluxes and producer
            B = T * obj.V;
            injector = max(obj.Qwells,0) .* dt/pv;  % injection flux * dt/pv
            
            ProductionSystem.Reservoir.State.S = ProductionSystem.Reservoir.State.S + (B * obj.f + injector);
        end
    end
end