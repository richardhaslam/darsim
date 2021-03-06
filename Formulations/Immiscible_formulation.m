% Immiscible Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_formulation < formulation
    properties
        Mobt
        Utot
        Qwells
        f
        df
        V
        Vr
        Ths
        UPc
        Ghs
    end
    methods
        %% constructor
        function obj = Immiscible_formulation()
            obj@formulation();
            obj.Tph = cell(2,1);
            obj.Gph = cell(2,1);
            obj.MatrixAssembler = matrix_assembler();
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
            %% 1. Reservoir Properties and Derivatives
            obj.drhodp = FluidModel.ComputeDrhoDp(ProductionSystem.Reservoir.State);
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.dMob = FluidModel.ComputeDMobDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.dPc = FluidModel.ComputeDPcDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            %% 2. Fractures Properties and Derivatives
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                obj.drhodp = [obj.drhodp; FluidModel.ComputeDrhoDp(ProductionSystem.FracturesNetwork.Fractures(f).State) ];
                obj.Mob = [obj.Mob; FluidModel.ComputePhaseMobilities(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('S_1').Value)];
                obj.dMob = [obj.dMob; FluidModel.ComputeDMobDS(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('S_1').Value)];
                obj.dPc = [obj.dPc; FluidModel.ComputeDPcDS(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('S_1').Value)];
            end
        end
        %% Methods for FIM Coupling
        function [Residual, RHS] = BuildMediumResidual(obj, Medium, Grid, dt, State0, Index, qw, qf, f, ph)
            % Create local variables
            N = Grid.N;
            pv = Medium.Por .* Grid.Volume;
            
            % Copy values in local variables
            s_old = State0.Properties(['S_', num2str(ph)]).Value(Index.Start:Index.End);
            rho_old = State0.Properties(['rho_', num2str(ph)]).Value(Index.Start:Index.End);
            P = Medium.State.Properties(['P_', num2str(ph)]).Value;
            s = Medium.State.Properties(['S_', num2str(ph)]).Value;
            rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            depth = Grid.Depth;
            
            % Accumulation Term
            if length(pv/dt) == 1
                AS = speye(N)*pv/dt;
            else
                AS = spdiags(pv/dt ,0,N,N);
            end
            
            % RESIDUAL
            RHS = qw(Index.Start:Index.End, ph);
            Residual = AS*(rho .* s - rho_old .* s_old)... % accummulation
                + obj.Tph{ph, 1+f} * P... % pressure
                - obj.Gph{ph, 1+f} * depth... % gravity
                - qw(Index.Start:Index.End, ph)... % wells
                - qf(Index.Start:Index.End, ph); % frac-matrix
        end
        function [Residual, RHS] = BuildFullResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute vector of qs
            qw = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            qf = zeros(DiscretizationModel.N, obj.NofPhases);
            if ProductionSystem.FracturesNetwork.Active
                qf = ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel);
            end
            % Transmissibility of reservoir
            for i=1:obj.NofPhases
                [obj.Tph{i, 1}, obj.Gph{i, 1}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{i, 1}, obj.Mob(1:DiscretizationModel.ReservoirGrid.N, i), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(i)]).Value, obj.GravityModel.RhoInt{i, 1});
            end
            
            % Initialise residual vector (Nph * N, 1)
            Nt = DiscretizationModel.N;
            Residual = zeros(Nt * obj.NofPhases, 1);
            RHS = zeros(Nt * obj.NofPhases, 1);
            Nm = DiscretizationModel.ReservoirGrid.N;
            
            for ph=1:obj.NofPhases
                % BuildResidual for Reservoir
                Index.Start = 1;
                Index.End = Nm;
                [Residualm, RHSm] = BuildMediumResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, qw, qf, 0, ph);
                Residual((ph-1)*Nt + Index.Start: (ph-1)*Nt + Index.End) = Residualm;
                RHS(     (ph-1)*Nt + Index.Start: (ph-1)*Nt + Index.End) = RHSm;
                
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                    % Transmissibility of fractures cells
                    [obj.Tph{ph, 1+f}, obj.Gph{ph, 1+f}] = ...
                        obj.MatrixAssembler.TransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph),...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, obj.GravityModel.RhoInt{ph, 1+f});
                    % BuildResidual for Fractures
                    Residual_frac_f = BuildMediumResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, qw, qf, f, ph);
                    Residual((ph-1)*Nt + Index.Start: (ph-1)*Nt + Index.End) = Residual_frac_f;
                end
            end
        end
        function [Jp, JS] = BuildMediumJacobian(obj, Medium, Wells, Grid, dt, Index, f, ph)
            % Create local variables
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            N_Face = length(obj.UpWind{ph,1+f});
            pv = Grid.Volume .* Medium.Por;
            
            rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            s = Medium.State.Properties(['S_', num2str(ph)]).Value;
            
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            % 1.a Pressure Block
            Jp = obj.Tph{ph,1+f};
            switch class(Grid)
                case('corner_point_grid')
                    nc = Grid.N;
                    nf = length(Grid.Trans);
                    C = [ Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index , Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index ];
                    D1 = [ -double(obj.U{ph,1+f}>=0)+double(obj.U{ph,1+f}<0)  , double(obj.U{ph,1+f}>=0)-double(obj.U{ph,1+f}<0)]; D1(D1==1)=0;
                    D2 = [ -double(obj.U{ph,1+f}>=0)+double(obj.U{ph,1+f}<0)  , double(obj.U{ph,1+f}>=0)-double(obj.U{ph,1+f}<0)];
                    UpwindPermutation1 = sparse([(1:nf)'; (1:nf)'], C, D1, nf, nc)';
                    UpwindPermutation2 = sparse([(1:nf)'; (1:nf)'], C, D2, nf, nc)';
                    
                    % 1.b: compressibility part
                    Mob_drhodp = obj.Mob(Index.Start:Index.End, ph) .* obj.drhodp(Index.Start:Index.End, ph);
                    dMobUpwind = obj.UpWind{ph,1+f} * Mob_drhodp;
                    FluxDerivative = dMobUpwind .* abs(obj.U{ph,1+f});
                    acc = pv/dt .* obj.drhodp(Index.Start:Index.End,ph) .* s;
                    Jp = Jp + spdiags(acc,0,N,N) + (UpwindPermutation1 * spdiags(FluxDerivative,0,N_Face,N_Face) * UpwindPermutation2')';
                    
                    % 2. Saturation Block
                    dMob_rho = obj.dMob(Index.Start:Index.End, ph) .* rho;
                    dMobUpwind = obj.UpWind{ph,1+f} * dMob_rho;
                    FluxDerivative = dMobUpwind .* abs(obj.U{ph,1+f});
                    v = (-1)^(ph+1) .* pv/dt .* rho;
                    JS = spdiags(v,0,N,N) + (UpwindPermutation1 * spdiags(FluxDerivative,0,N_Face,N_Face) * UpwindPermutation2')';
                    
                case('cartesian_grid')
                    % 1.b: compressibility part
                    dMupx = obj.UpWind{ph,1+f}.x * ( obj.Mob(Index.Start:Index.End, ph) .* obj.drhodp(Index.Start:Index.End, ph) );
                    dMupy = obj.UpWind{ph,1+f}.y * ( obj.Mob(Index.Start:Index.End, ph) .* obj.drhodp(Index.Start:Index.End, ph) );
                    dMupz = obj.UpWind{ph,1+f}.z * ( obj.Mob(Index.Start:Index.End, ph) .* obj.drhodp(Index.Start:Index.End, ph) );
                    
                    X1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                    X2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                    Y1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                    Y2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                    Z1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                    Z2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
                    acc = pv/dt .* obj.drhodp(Index.Start:Index.End,ph) .* s;
                    
                    DiagVecs = [-Z2, -Y2, -X2, Z2+Y2+X2-Z1-Y1-X1+acc, X1, Y1, Z1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    Jp = Jp + spdiags(DiagVecs, DiagIndx, N, N);
                    
                    % 2. Saturation Block
                    dMupx = obj.UpWind{ph,1+f}.x * ( obj.dMob(Index.Start:Index.End, ph) .* rho );
                    dMupy = obj.UpWind{ph,1+f}.y * ( obj.dMob(Index.Start:Index.End, ph) .* rho );
                    dMupz = obj.UpWind{ph,1+f}.z * ( obj.dMob(Index.Start:Index.End, ph) .* rho );
                    
                    % Construct JS block
                    X1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                    X2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                    Y1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                    Y2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                    Z1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                    Z2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
                    v = (-1)^(ph+1) .* pv/dt .* rho;
                    DiagVecs = [-Z2, -Y2, -X2, Z2+Y2+X2-Z1-Y1-X1+v, X1, Y1, Z1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    JS = spdiags(DiagVecs,DiagIndx,N,N);
            end
            
            %Add capillarity
            if ph == 1
                JS = JS - Jp * spdiags(obj.dPc, 0, N, N);
            end
            
            % Add Wells
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [Jp, JS] = obj.AddWellsToJacobian(Jp, JS, Medium.State, Wells, Medium.K(:,1), ph);
            end
        end
        function Jacobian = BuildFullJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %% Jacobian's assembly
            %      pm     pf   |   sm      sf
            % m| J1m_pm J1m_pf | J1m_sm J1m_sm | | dpm | = | R1m |
            % f| J1f_pm J1f_pf | J1f_sm J1f_sf | | dpf |   | R1f |
            %  |---------------|---------------|
            % m| J2m_pm J2m_pf | J2m_sm J2m_sf | | dsm |   | R2m |
            % f| J2f_pm J2f_pf | J2f_sm J2f_sf | | dsf |   | R2f |
            %  |---------------|---------------|
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            Reservoir = ProductionSystem.Reservoir;
            Fractures = ProductionSystem.FracturesNetwork.Fractures;
            Wells = ProductionSystem.Wells;
            % Global variables
            if ProductionSystem.FracturesNetwork.Active
                FineGrid = [DiscretizationModel.ReservoirGrid; DiscretizationModel.FracturesGrid.Grids];
            else
                FineGrid = DiscretizationModel.ReservoirGrid;
            end
            P = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'P_'); % useful for cross connections assembly
            rho = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'rho_'); % useful for cross connections assembly
            
            Jph = cell(obj.NofPhases, 1);
            for ph=1:obj.NofPhases
                Jph{ph} = sparse(Nt, obj.NofPhases*Nt);
                %% Jacobian of the reservoir
                Index.Start = 1;
                Index.End = Nm;

                [Jp_res, JS_res] = BuildMediumJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph);
                
                %% Jacobian of the fractures
                Jp = [];  JS = [];
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global( Nm, f, 1     );
                    Index.End   = DiscretizationModel.Index_Local_to_Global( Nm, f, Nf(f) );
                    [Jp_frac, JS_frac] = BuildMediumJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    Jp  = blkdiag(Jp, Jp_frac);
                    JS  = blkdiag(JS, JS_frac);
                end
                
                Jp  = blkdiag(Jp_res,Jp);
                JS  = blkdiag(JS_res,JS);

                %% ADD frac-matrix and frac-frac connections
                for c = 1:length(DiscretizationModel.CrossConnections)
                    if isempty(DiscretizationModel.CrossConnections(c).Cells),  continue;  end
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    %% 1. Pressure block
                    Jp_conn = - T_Geo .* (...
                             UpWind(:, ph) .* obj.Mob(j, ph) .* (rho(j,ph) + obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) )) + ...
                            ~UpWind(:, ph) .* obj.Mob(i, ph) .* (rho(i,ph) + obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) )) );
                    % frac - mat or frac1 - frac2
                    Jp(i, j) = Jp_conn;
                    Jp(i, i) = Jp(i, i) - sum(Jp_conn);
                    % mat-frac or frac2 - frac1
                    Jp(j, i) = Jp_conn';
                    Jp(sub2ind([Nt, Nt], j, j)) = Jp(sub2ind([Nt, Nt], j, j)) - Jp_conn;
                    %% 2. Saturation block
                    JS1_conn = T_Geo .*  UpWind(:, ph) .* (P(i, ph) - P(j, ph)).* (rho(j, ph) .* obj.dMob(j, ph));
                    JS2_conn = T_Geo .* ~UpWind(:, ph) .* (P(i, ph) - P(j, ph)).* (rho(i, ph) .* obj.dMob(i, ph));
                    % frac - mat or frac1 - frac2
                    JS(i, j) = JS1_conn;
                    JS(i, i) = JS(i, i) + sum(JS2_conn);
                    % mat-frac or frac2 - frac1
                    JS(j, i) = -JS2_conn;
                    % diag of mat or frac2
                    JS(sub2ind([Nt, Nt], j, j)) = JS(sub2ind([Nt, Nt], j, j)) - JS1_conn;
                end
                
                % Combine the Jacobian blocks
                Jph{ph} = horzcat(Jp, JS);
            end
            Jacobian = Jph{1};
            %Jacobian_noconn = Jphnoconn{1};
            for i=2:obj.NofPhases
                Jacobian = vertcat(Jacobian, Jph{i});
                %Jacobian_noconn = vertcat(Jacobian_noconn, Jphnoconn{i});
            end
        end
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                Nm = DiscretizationModel.ReservoirGrid.N;
                Nt = DiscretizationModel.N;
                deltaP = delta(1:Nt);
                deltaS = delta(Nt+1: end);
                %% 1. Update matrix
                % Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(obj.NofPhases)]);
                Pm.update(deltaP(1:Nm));
                DeltaLast = zeros(Nm, 1);
                for ph = 1:obj.NofPhases-1
                    Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(ph)]);
                    DeltaS = deltaS(1:Nm);
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
                    EP = Nm;
                    Nf = DiscretizationModel.FracturesGrid.N;
                    for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                        IP = EP+1;
                        EP = IP + Nf(f) - 1;
                        % Update Pressure
                        Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(obj.NofPhases)]);
                        Pf.update(deltaP(IP:EP));
                        DeltaLast = zeros(Nf(f), 1);
                        for ph = 1:obj.NofPhases-1
                            IS = Nt*(ph-1) + IP;
                            ES = Nt*(ph-1) + EP;
                            Sf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['S_', num2str(ph)]);
                            Sf.update(deltaS(IS:ES));
                            Sf.Value = max(Sf.Value, 0);
                            Sf.Value = min(Sf.Value, 1);
                            DeltaLast = DeltaLast + deltaS(IS:ES);
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
        function qw = ComputeSourceTerms(obj, N, Wells)
            qw = zeros(N, obj.NofPhases);    
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj(i).Cells;
                qw(c, :) = Wells.Inj(i).QPhases(:,:);
            end
            %Producers
            for i=1:Wells.NofProd
                c = Wells.Prod(i).Cells;
                qw(c, :) = Wells.Prod(i).QPhases(:,:);
            end
        end
        function qf = ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel)
            qf = zeros(DiscretizationModel.N, obj.NofPhases);
            Nm = DiscretizationModel.ReservoirGrid.N;
            % Global variables
            P = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid; DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'P_'); % useful for cross connections assembly
            rho = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid; DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'rho_'); % useful for cross connections assembly
            for ph=1:obj.NofPhases
                % fill in qf
                for c=1:length(DiscretizationModel.CrossConnections)
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    i = c + Nm;
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    qf(i, ph) = qf(i, ph) + sum(T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* (P(j, ph) - P(i, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* (P(j, ph) - P(i, ph))) );
                    qf(j, ph) = qf(j, ph) + T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* (P(i, ph) - P(j, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* (P(i, ph) - P(j, ph)) );
                end
            end
        end
        function [Jp, JS] = AddWellsToJacobian(obj, Jp, JS, State, Wells, K, ph)
            % Define Local handles
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %p = State.Properties('P_2').Value;
            %rho = State.Properties(['rho_', num2str(ph)]).Value;
          
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                dQdp = Inj(i).ComputeWellMassFluxDerivativeWithRespectToPressure(K, obj.NofPhases);
                for j=1:length(a)
                    Jp(a(j),a(j)) = Jp(a(j),a(j)) - dQdp(j, ph);
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                dQdp = Prod(i).ComputeWellMassFluxDerivativeWithRespectToPressure(State, K, obj.Mob, obj.drhodp, zeros(size(K,1),2), obj.NofPhases);
                dQdS = Prod(i).ComputeWellMassFluxDerivativeWithRespectToSaturation(State, K, obj.dMob, obj.NofPhases);
                for j=1:length(b)
                    Jp(b(j),b(j)) = Jp(b(j),b(j)) - dQdp(j, ph);                    
                    JS(b(j),b(j)) = JS(b(j),b(j)) - dQdS(j, ph);
                end
            end
        end
        function AverageMassOnCoarseBlocks(obj, ProductionSystem, FineGrid, FluidModel, R)
            S = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'S_');
            
            % Perform Average for ADM
            delta = zeros(sum([FineGrid.N]), obj.NofPhases);
            for ph = 1:obj.NofPhases - 1
                S_rest = R * S(:,ph);
                Sav = R' * (S_rest ./ sum(R, 2));
                delta(:, ph) = Sav - S(:,ph);
            end
            delta(:, end) = -sum(delta(:,1:end-1), 2);
            
            for ph=1:obj.NofPhases
                Start=1;
                End = FineGrid(1).N;
                S = ProductionSystem.Reservoir.State.Properties(['S_', num2str(ph)]);
                S.update(delta(Start:End, ph));
                for frac = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    Start = End + 1;
                    End = Start + FineGrid(frac+1).N - 1;
                    S = ProductionSystem.FracturesNetwork.Fractures(frac).State.Properties(['S_', num2str(ph)]);
                    S.update(delta(Start:End, ph));
                end
            end
            % Here FluidModel is useless but for compositional formulations it
            % is important
        end
        function CFL = ComputeCFLNumber(obj, ProductionSystem, DiscretizationModel, dt)
            N = DiscretizationModel.ReservoirGrid.N;      
            pv = ProductionSystem.Reservoir.Por.*DiscretizationModel.ReservoirGrid.Volume;
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
                [obj.Tph{i}, obj.Gph{i}] = obj.MatrixAssembler.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, obj.UpWind{i, 1}, obj.Mob(1:N,i), rho(:,i), obj.GravityModel.RhoInt{i, 1});
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
                Mass(:,i) = rho(:,i) .* S(:,i) .* pv;
            end
            Mass = max(Mass, 1e-4);
            %ThroughPut(ThroughPut < 1e-4) = 0;
            Ratio = ThroughPut ./ Mass;
            CFL = dt * max(max(Ratio));
        end
        %% Methods for Sequential Coupling
        function ComputeTotalMobility(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.Mobt = sum(obj.Mob, 2);
        end
        function UpdateFractionalFlow(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.f = obj.Mob(:,1) ./ (obj.Mob(:,1) + obj.Mob(:, 2));
        end
        function UpdateCapillaryPressure(obj, ProductionSystem, FluidModel)
            FluidModel.ComputePc(ProductionSystem.Reservoir.State);
            obj.dPc = FluidModel.ComputeDPcDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
        end
        function ComputeDfDS(obj, ProductionSystem, FluidModel)
            obj.dMob = FluidModel.ComputeDMobDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            num = obj.Mob(:, 1);
            dnum = obj.dMob(:,1);
            den = sum(obj.Mob, 2);
            dden = sum(obj.dMob, 2);
            obj.df = (dnum .* den - dden .* num) ./ den.^2;
        end
        function dfdSdS = ComputeDfDSDS(obj, s, rho, FluidModel)
           % f = (rho1 * Mob1) / (rho1 * Mob1 + rho2 * Mob2);
           % df = (rho1 * dMob1 * (rho1*Mob1 + rho2*Mob2) - (rho1*dMob1 + rho2*dMob2) * rho1*Mob1) / (rho1 * Mob1 + rho2 * Mob2)^2
           % Num = A - B;
           % A = (rho1 * dMob1 * (rho1*Mob1 + rho2*Mob2)
           % dA = rho1 * ddMob1 * (rho1*Mob1 + rho2*Mob2) + rho1 * dMob1 * (rho1*dMob1 + rho2 * dMob2)
           % B = (rho1*dMob1 + rho2*dMob2) * rho1*Mob1
           % dB = (rho1*ddMob1 + rho2 * ddMob2) *  rho1*Mob1 + (rho1*dMob1 + rho2*dMob2) * rho1*dMob1 
           % Den = (rho1 * Mob1 + rho2 * Mob2)^2;
           Mob = FluidModel.ComputePhaseMobilities(s);
           dMob = FluidModel.ComputeDMobDS(s);
           ddMob = FluidModel.ComputeDMobDSDS(s);
         
           A = rho(:, 1) .* dMob(:, 1) .* sum(rho.*Mob, 2);
           dA = rho(:, 1) .* ddMob(:,1) .* sum(rho.*Mob, 2) + rho(:,1) .* dMob(:,1) .*  sum(rho.*dMob, 2);
           B = sum(rho.*dMob, 2) .* rho(:, 1) .* Mob(:,1);
           dB = sum(rho.*ddMob, 2) .*  rho(:, 1).*Mob(:, 1) + sum(rho.*dMob, 2) .* rho(:, 1) .* dMob(:, 1);
           num = A - B;
           dnum = dA - dB;
           den = sum(rho.*Mob, 2).^2;
           dden = sum(rho.*Mob, 2) .* sum(rho.*dMob, 2);  
           
           dfdSdS = (dnum .* den - dden .* num) ./ den.^2;
           
        end
        function Residual = BuildPressureResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            
            %%%% fix wells
            % Compute vector of qs
            qw = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            qf = zeros(DiscretizationModel.N, obj.NofPhases);
            if ProductionSystem.FracturesNetwork.Active
                qf = ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel);
            end
            % Transmissibility of reservoir
            for i=1:obj.NofPhases
                [obj.Tph{i, 1}, obj.Gph{i, 1}] = ...
                    obj.MatrixAssembler.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, obj.UpWind{i, 1}, obj.Mob(1:DiscretizationModel.ReservoirGrid.N, i), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(i)]).Value, obj.GravityModel.RhoInt{i, 1});
            end
            
            % Initialise residual vector (Nph * N, 1)
            Nt = DiscretizationModel.N;
            Nm = DiscretizationModel.ReservoirGrid.N;
            Residual = zeros(Nt, 1);
            
            % BuildResidual for Reservoir
            Index.Start = 1;
            Index.End = Nm;
            phi = ProductionSystem.Reservoir.Por;

            Residual(Index.Start:Index.End) = MediumPressureResidual(obj, DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir.State, State0, dt, phi, qw, qf, Index, 0);
            % Fractures
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + Nf(f) - 1;
                    % Transmissibility of fractures cells
                    for i=1:obj.NofPhases
                        [obj.Tph{i, 1+f}, obj.Gph{i, 1+f}] = ...
                            obj.MatrixAssembler.TransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), obj.UpWind{i, 1+f}, obj.Mob(Index.Start:Index.End, i),...
                            ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(i)]).Value, obj.GravityModel.RhoInt{i, f+1});
                    end
                    % BuildResidual for Fractures
                    
                    phi = ProductionSystem.FracturesNetwork.Fractures(f).Por;
                    Residual(Index.Start:Index.End) = MediumPressureResidual(obj, DiscretizationModel.FracturesGrid.Grids(f), ProductionSystem.FracturesNetwork.Fractures(f).State, State0, dt, phi, qw, qf, Index, f);
                end
            end
        end
        function Residual = MediumPressureResidual(obj, Grid, State, State0, dt, phi, qw, qf, Index, f)
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
                obj.Gph{i, f+1} = bsxfun(@rdivide, obj.Gph{i, f+1}, rho(:, i));
            end 
            depth = Grid.Depth;
            Residual = zeros(N, 1);
            
            for i = 1:obj.NofPhases
                Residual(:) = ...
                    Residual(:) - ...
                    pv/dt * (rho_old(:, i) .* s_old(:, i) ./ rho(:, i)) ...
                    + obj.Tph{i, f+1} * P(:, i) ...
                    - obj.Gph{i, f+1} * depth ...
                    - qw(Index.Start:Index.End, i) ./ rho(:, i)...
                    - qf(Index.Start:Index.End, i) ./ rho(:, i);
            end
            Residual = Residual + pv/dt;
        end
        function A = BuildPressureMatrix(obj, ProductionSystem, DiscretizationModel, dt, State0)
            %% Pressure Matrix's assembly
            % m|   Pmm   |   Pmf   |
            %  |---------|---------|
            % f|   Pfm   |   Pff   |
            
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            Nm = DiscretizationModel.ReservoirGrid.N;
            % Global variables
            if ProductionSystem.FracturesNetwork.Active
                FineGrid = [DiscretizationModel.ReservoirGrid; DiscretizationModel.FracturesGrid.Grids];
            else
                FineGrid = DiscretizationModel.ReservoirGrid;
            end
            % P = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'P_'); % useful for cross connections assembly
            rho = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'rho_'); % useful for cross connections assembly        
            %% 1. Reservoir
            Index.Start = 1;
            Index.End = Nm;
            A = obj.Build_Medium_PressureMatrix(ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, 0);
            % Add wells
            W = obj.AddWellsToPressureSystem(DiscretizationModel.ReservoirGrid.N, ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K(:,1), rho(1:Nm,:));
            A = A + W;
            %% 2. Fractures
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                    Index.End   = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, Nf(f));
                    Af = obj.Build_Medium_PressureMatrix(ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, f);
                    A = blkdiag(A, Af);
                end
            end
            %% 3.  Add matrix-frac and frac-frac connections
            for ph=1:obj.NofPhases
                for c = 1:length(DiscretizationModel.CrossConnections)
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    A(i,j) = A(i,j) - T_Geo' .* obj.Mob(j,ph)' .* rho(j,ph)';
                    A(i,i) = A(i,i) + sum(T_Geo.* obj.Mob(j,ph).* rho(j,ph));
                    A(j,i) = A(j,i) - T_Geo.* obj.Mob(i,ph).* rho(i,ph);
                    A(sub2ind(size(A), j, j)) = A(sub2ind(size(A), j, j)) + T_Geo.* obj.Mob(i,ph).* rho(i,ph);
                end
            end
        end
        function J = Build_Medium_PressureMatrix(obj, Medium, Grid, dt, State0, Index, f)
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
                s_old(:, i) = State0.Properties(['S_', num2str(i)]).Value(Index.Start:Index.End);
                rho_old(:, i) = State0.Properties(['rho_', num2str(i)]).Value(Index.Start:Index.End);
                rho(:, i) = Medium.State.Properties(['rho_', num2str(i)]).Value;
            end 
            
            A = sparse(N,N);
            C = sparse(N,N);
            for i=1:obj.NofPhases
                %% 1. Accummulation term
                vec = - pv/dt * (- rho_old(:, i) .* s_old(:, i) .* obj.drhodp(Index.Start:Index.End,i) ./ rho(:, i).^2) ;
                % Alternative
                %vec = pv/dt * obj.drhodp(Start:End,i);
                C = C + spdiags(vec, 0, N, N);
                
                %% 2. Convective term 
                % a. transmissibility matrix
                A = A + obj.Tph{i, f+1};
                % b. Derivative of transmissibility with respect to p
                dMupx = obj.UpWind{i, f+1}.x*(obj.Mob(Index.Start:Index.End, i) .* obj.drhodp(Index.Start:Index.End,i));
                dMupy = obj.UpWind{i, f+1}.y*(obj.Mob(Index.Start:Index.End, i) .* obj.drhodp(Index.Start:Index.End,i));
                dMupz = obj.UpWind{i, f+1}.z*(obj.Mob(Index.Start:Index.End, i) .* obj.drhodp(Index.Start:Index.End,i));
                
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
                [dq(c, :), ~] = Inj(i).dQPhasesdPdS(K, obj.NofPhases);
            end
            
            %Producers
            for i=1:length(Prod)
                c = Prod(i).Cells;
                [dq(c, :), ~] = Prod(i).dQPhasesdPdS(State, K, obj.Mob, obj.dMob, obj.drhodp, obj.NofPhases);
            end
            q = sparse(obj.ComputeSourceTerms(N, Wells));
            dqt = sum((dq .* rho - obj.drhodp(1:N,:).* q)./rho.^2, 2);
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
            N  = DiscretizationModel.ReservoirGrid.N;

            rho = zeros(N, obj.NofPhases);
            for i=1:obj.NofPhases
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end 
            
            % Wells total fluxes
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            obj.Qwells = sum(q./rho, 2);
            
            switch class(DiscretizationModel.ReservoirGrid)
                case('corner_point_grid')
                    N_Faces = length(obj.U{1,1});
                    obj.Utot = zeros(N_Faces,1);
                    for ph = 1:obj.NofPhases
                        MobUpwind = obj.UpWind{ph,1} * obj.Mob(1:N, ph);
                        obj.Utot = obj.Utot + obj.U{ph,1} .* MobUpwind;
                    end
                    if ProductionSystem.Reservoir.State.Properties('V_tot').Plot
                       % Not yet implemented
                    end
                case('cartesian_grid')
                    obj.Utot.x = zeros(Nx+1, Ny, Nz);
                    obj.Utot.y = zeros(Nx, Ny+1, Nz);
                    obj.Utot.z = zeros(Nx, Ny, Nz+1);
                    for ph = 1:obj.NofPhases
                        obj.Utot.x(2:Nx+1,:,:) = obj.Utot.x(2:Nx+1,:,:) + obj.U{ph, 1}.x(2:Nx+1,:,:) .* reshape(obj.UpWind{ph,1}.x *  (obj.Mob(1:N, ph)), Nx, Ny, Nz);  %- Ucap.x(2:Nx,:);
                        obj.Utot.y(:,2:Ny+1,:) = obj.Utot.y(:,2:Ny+1,:) + obj.U{ph, 1}.y(:,2:Ny+1,:) .* reshape(obj.UpWind{ph,1}.y *  (obj.Mob(1:N, ph)), Nx, Ny, Nz);  %- Ucap.y(:,2:Ny);
                        obj.Utot.z(:,:,2:Nz+1) = obj.Utot.z(:,:,2:Nz+1) + obj.U{ph, 1}.z(:,:,2:Nz+1) .* reshape(obj.UpWind{ph,1}.z *  (obj.Mob(1:N, ph)), Nx, Ny, Nz);  %- Ucap.y(:,2:Ny);
                    end
                    if ProductionSystem.Reservoir.State.Properties('V_tot').Plot
                        % Compute average between the 2 interfaces
                        ux = (obj.Utot.x(1:Nx,:,:) + obj.Utot.x(2:Nx+1,:,:))/2;
                        uy = (obj.Utot.y(:,1:Ny,:) + obj.Utot.y(:,2:Ny+1,:))/2;
                        uz = (obj.Utot.z(:,:,1:Nz) + obj.Utot.z(:,:,2:Nz+1))/2;
                        % Copy the values in the right objects
                        Ucentres = ProductionSystem.Reservoir.State.Properties('V_tot');
                        Ucentres.Value(:,1) = reshape(ux, N, 1);
                        Ucentres.Value(:,2) = reshape(uy, N, 1);
                        Ucentres.Value(:,3) = reshape(uz, N, 1);
                    end
            end
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
            ux = reshape(obj.Utot.x(1:Nx,:,:) - obj.Utot.x(2:Nx+1,:,:), Grid.N, 1);
            uy = reshape(obj.Utot.y(:,1:Ny,:) - obj.Utot.y(:,2:Ny+1,:), Grid.N, 1);
            uz = reshape(obj.Utot.z(:,:,1:Nz) - obj.Utot.z(:,:,2:Nz+1), Grid.N, 1);
            Balance = ux + uy + uz + obj.Qwells;
            if norm(Balance/maxU, inf) > 1e-5
                conservative = 0;
            end
        end
        function Residual = BuildTransportResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            
            N = DiscretizationModel.ReservoirGrid.N;
           
            % Initialise local objects
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            s = ProductionSystem.Reservoir.State.Properties('S_1').Value;
            s_old = State0.Properties('S_1').Value;      
            rho = ProductionSystem.Reservoir.State.Properties('rho_1').Value;
            % viscous fluxes matrix
            obj.V = obj.MatrixAssembler.ViscousMatrix(DiscretizationModel.ReservoirGrid, obj.Qwells, obj.Utot); 
            vect = obj.Mob(1:N, 2) .* obj.f;
            
            Pc = ProductionSystem.Reservoir.State.Properties('Pc').Value;
            [obj.Ths, obj.UPc, obj.Ghs] =  obj.MatrixAssembler.TransmissibilityforS(DiscretizationModel.ReservoirGrid, obj.UpWind{1, 1}, vect, rho, Pc, obj.GravityModel.RhoInt);
            
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            % Compute residual
            Residual = pv/dt .* rho .* (s - s_old)  - max(q(:, 1), 0) - obj.V * (obj.f .* rho) - ...
                       obj.Ths * (Pc)- obj.Ghs * (depth);
        end
        function Jacobian = BuildTransportJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            % Build Transport Jacobian
            N = DiscretizationModel.ReservoirGrid.N;
            
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            rho = ProductionSystem.Reservoir.State.Properties('rho_1').Value;
            
            % 1.a: saturation partial derivative in time + total velocity
            
            D = spdiags(pv/dt*rho, 0, N, N);
            Jacobian = D - obj.V * spdiags(obj.df.*rho,0,N,N);
            Jacobian = Jacobian - obj.Ths * spdiags(obj.dPc,0,N,N);
       
            % 1.a: capillarity part -> diffusive term
            dMupx = obj.UpWind{1,1}.x*(obj.Mob(1:N, 2) .* obj.df .* rho);
            dMupy = obj.UpWind{1,1}.y*(obj.Mob(1:N, 2) .* obj.df .* rho);
            dMupz = obj.UpWind{1,1}.z*(obj.Mob(1:N, 2) .* obj.df .* rho);
                      
            
            vecX1 = min(reshape(obj.UPc.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
            vecX2 = max(reshape(obj.UPc.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
            vecY1 = min(reshape(obj.UPc.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
            vecY2 = max(reshape(obj.UPc.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
            vecZ1 = min(reshape(obj.UPc.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
            vecZ2 = max(reshape(obj.UPc.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
           
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1, vecX1, vecY1, vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            Jacobian  = Jacobian - spdiags(DiagVecs, DiagIndx, N, N);
            
            dMupx = obj.UpWind{1,1}.x*(obj.dMob(1:N, 2) .* obj.f .* rho);
            dMupy = obj.UpWind{1,1}.y*(obj.dMob(1:N, 2) .* obj.f .* rho);
            dMupz = obj.UpWind{1,1}.z*(obj.dMob(1:N, 2) .* obj.f .* rho);
 
            vecX1 = min(reshape(obj.UPc.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
            vecX2 = max(reshape(obj.UPc.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
            vecY1 = min(reshape(obj.UPc.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
            vecY2 = max(reshape(obj.UPc.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
            vecZ1 = min(reshape(obj.UPc.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
            vecZ2 = max(reshape(obj.UPc.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
           
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1, vecX1, vecY1, vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            Jacobian  = Jacobian - spdiags(DiagVecs, DiagIndx, N, N);
            
        end
        function delta = UpdateSaturation(obj, ProductionSystem, delta, FluidModel, DiscretizationModel)
            s_old = ProductionSystem.Reservoir.State.Properties('S_1').Value;
            rho = 0 * obj.Mob;
            for i=1:obj.NofPhases
                rho(:,i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end
            
            % Update
            snew = s_old + delta;
            
            % Remove values that are not physical
            snew = max(snew, 0);
            snew = min(snew, 1);
            
            % FLUX CORRECTION - PATRICK
            Ddf_old = obj.ComputeDfDSDS(s_old, rho, FluidModel);
            Ddf = obj.ComputeDfDSDS(snew, rho, FluidModel);           
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
        function CFL = ComputeCFLNumberTransport(obj, DiscretizationModel, ProductionSystem, dt)
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            dx = DiscretizationModel.ReservoirGrid.dx;
            dy = DiscretizationModel.ReservoirGrid.dy;
            dz = DiscretizationModel.ReservoirGrid.dz;
            maxdf = max(abs(obj.df));
            maxUx = max(max(max(abs(obj.Utot.x))));
            maxUy = max(max(max(abs(obj.Utot.y))));
            maxUz = max(max(max(abs(obj.Utot.z))));
            CFL = dt * maxdf * (maxUx + maxUy + maxUz) / pv;
        end
        %% Explicit transport solver 
        function delta = UpdateSaturationExplicitly(obj, ProductionSystem, DiscretizationModel, dt)
            % 0. Initialise
            pv = ProductionSystem.Reservoir.Por .* DiscretizationModel.ReservoirGrid.Volume;
            N = DiscretizationModel.ReservoirGrid.N;
            
            obj.V = obj.MatrixAssembler.ViscousMatrix(DiscretizationModel.ReservoirGrid, obj.Qwells, obj.Utot); 

            % 1. Solve
            T = spdiags(dt/pv*ones(N,1),0,N,N);    % dt/pv * Cell Fluxes and producer
            B = T * obj.V;
            injector = max(obj.Qwells, 0) .* dt/pv;  % injection flux * dt/pv
            
            
            S = ProductionSystem.Reservoir.State.Properties('S_1');
            s_old = S.Value;
            vect =  obj.Mob(1:N, 2) .* obj.f; 
            Pc = ProductionSystem.Reservoir.State.Properties('Pc').Value; 
            
            [obj.Ths, ~] = obj.MatrixAssembler.TransmissibilityforS(DiscretizationModel.ReservoirGrid, obj.UpWind{1, 1}, vect, ones(size(s_old)), Pc);
            Pcap = T * obj.Ths;
            
            
            delta = S.Value - s_old + (B * obj.f  + injector) + Pcap * (Pc);
            
            % Now update values
            Nm = DiscretizationModel.ReservoirGrid.N;
            DeltaLast = zeros(Nm, 1);
            for ph = 1:obj.NofPhases-1
                Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(ph)]);
                DeltaS = delta(1:Nm);
                Sm.update(DeltaS);
                % Remove values that are not physical
                DeltaLast = DeltaLast + DeltaS;
            end
            Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(obj.NofPhases)]);
            Sm.update(-DeltaLast);

            if max(Sm.Value) > 1.01 || min(Sm.Value) < -0.01
                error('DARSim2 error: Unphysical Saturation value: try reducing CFL number!')
            end
        end
    end
end