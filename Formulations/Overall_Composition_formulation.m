% Overall Composition Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 7 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Overall_Composition_formulation < Compositional_formulation
    properties
        dxdz
        dxdp
        drhoTdp
        drhodz
        drhoTdz
        dMobdp
        dMobdz
    end
    methods
        function obj = Overall_Composition_formulation(n_components)
            obj@Compositional_formulation(n_components);
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            % This is the bitchy part!! 
            [obj.dxdp, obj.dxdz] = FluidModel.DxDpDz(ProductionSystem.Reservoir.State, obj.SinglePhase);
            % %%%%%%%%%%%%%%%%%%%%%%%%
            obj.drhodp = FluidModel.DrhoDp(ProductionSystem.Reservoir.State, obj.SinglePhase);
            obj.drhodz = FluidModel.DrhoDz(ProductionSystem.Reservoir.State, obj.SinglePhase);
            dSdp = FluidModel.DSDp(ProductionSystem.Reservoir.State, obj.drhodp, -obj.dxdp(:,5));
            dSdz = FluidModel.DSDz(ProductionSystem.Reservoir.State, -obj.dxdz(:, end));
            obj.drhoTdz = FluidModel.DrhotDz(ProductionSystem.Reservoir.State, obj.drhodz, dSdz);
            obj.drhoTdp = FluidModel.DrhotDp(ProductionSystem.Reservoir.State,obj.drhodp, dSdp);
            obj.dMobdp = FluidModel.DMobDp(ProductionSystem.Reservoir.State, [dSdp, dSdp]); % I am not using it !
            obj.dMob = FluidModel.DMobDz(ProductionSystem.Reservoir.State, dSdz);
            obj.dPc = FluidModel.DPcDz(ProductionSystem.Reservoir.State, dSdz);
         end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            %Create local variables
            N = DiscretizationModel.ReservoirGrid.N;          
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            z_old = zeros(N, obj.NofComponents);
            rho_old = zeros(N, obj.NofPhases);
            P = zeros(N, obj.NofPhases);
            z = zeros(N, obj.NofComponents);
            rho = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofComponents*obj.NofPhases);
            
            % Copy values in local variables
            for i=1:obj.NofPhases
                rho_old(:,i) = State0.Properties(['rho_', num2str(i)]).Value;
                P(:, i) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end
            for i=1:obj.NofComponents
                z_old(:,i) = State0.Properties(['z_', num2str(i)]).Value;
                z(:, i) = ProductionSystem.Reservoir.State.Properties(['z_', num2str(i)]).Value;
                for j=1:obj.NofPhases
                    x(:,(i-1)*obj.NofPhases + j) = ProductionSystem.Reservoir.State.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            rhoT_old = State0.Properties('rhoT').Value;
            rhoT =  ProductionSystem.Reservoir.State.Properties('rhoT').Value;
            
            % Depths
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            %Accumulation term
            A = speye(N)*pv/dt;
            
            %Source terms
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            %%  Build Residual
            Residual = zeros(N*obj.NofComponents, 1);
            for i=1:obj.NofComponents
                % Total moles
                n = rhoT .* z(:,i);
                n_old = rhoT_old .* z_old(:,i);
                % Phase Transmissibilities
                obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, rho, obj.GravityModel.RhoInt, x(:,(i-1)*2+1:(i-1)*2+2), i); 
                Residual((i-1)*N+1:i*N) = A * (n - n_old) ...             % Accumulation term
                           + obj.Tph{i, 1} *  P(:,1) ...    % Convective term                
                           + obj.Tph{i, 2} *  P(:,2)...
                           + obj.Gph{i,1} * depth...       % Gravity
                           + obj.Gph{i,2} * depth...
                           - q(:,i);                       % Source/Sink  
            end            
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            
            % Initialise local variables
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N = DiscretizationModel.ReservoirGrid.N;
            %
            pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            z = zeros(N, obj.NofComponents);
            rho = zeros(N, obj.NofPhases);
            P = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofComponents*obj.NofPhases);
            
            % Copy values in local variables
            for i=1:obj.NofPhases
                P(:, i) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end
            for i=1:obj.NofComponents
                z(:, i) = ProductionSystem.Reservoir.State.Properties(['z_', num2str(i)]).Value;
                for j=1:obj.NofPhases
                    x(:,(i-1)*obj.NofPhases + j) = ProductionSystem.Reservoir.State.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            rhoT =  ProductionSystem.Reservoir.State.Properties('rhoT').Value;
            
            % Fill in block by block
            Jp = cell(obj.NofComponents, 1);
            Jz = cell(obj.NofComponents, obj.NofComponents-1);
            for i=1:obj.NofComponents
                %% 1. Component i pressure block
                % 1.a: divergence
                Jp{i} = obj.Tph{i,1}  + obj.Tph{i, 2};
                % 1.b: compressibility part
                dMupxPh1 = obj.UpWind(1).x * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodp(:,1) + obj.Mob(:, 1) .* obj.dxdp(:,(i-1)*2+1) .* rho(:,1));% + obj.dMobdp(:, 1) .* x(:,(i-1)*2+1) .* rho(:,1));
                dMupyPh1 = obj.UpWind(1).y * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodp(:,1) + obj.Mob(:, 1) .* obj.dxdp(:,(i-1)*2+1) .* rho(:,1));% + obj.dMobdp(:, 1) .* x(:,(i-1)*2+1) .* rho(:,1));
                dMupzPh1 = obj.UpWind(1).z * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodp(:,1) + obj.Mob(:, 1) .* obj.dxdp(:,(i-1)*2+1) .* rho(:,1));% + obj.dMobdp(:, 1) .* x(:,(i-1)*2+1) .* rho(:,1));
                dMupxPh2 = obj.UpWind(2).x * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drhodp(:,2) + obj.Mob(:, 2) .* obj.dxdp(:,(i-1)*2+2) .* rho(:,2));% + obj.dMobdp(:, 2) .* x(:,(i-1)*2+2) .* rho(:,2));
                dMupyPh2 = obj.UpWind(2).y * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drhodp(:,2) + obj.Mob(:, 2) .* obj.dxdp(:,(i-1)*2+2) .* rho(:,2));% + obj.dMobdp(:, 2) .* x(:,(i-1)*2+2) .* rho(:,2));
                dMupzPh2 = obj.UpWind(2).z * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drhodp(:,2) + obj.Mob(:, 2) .* obj.dxdp(:,(i-1)*2+2) .* rho(:,2));% + obj.dMobdp(:, 2) .* x(:,(i-1)*2+2) .* rho(:,2));
                
                vecX1 = min(reshape(obj.U(1).x(1:Nx,:,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U(2).x(1:Nx,:,:),N,1), 0).*dMupxPh2;
                vecX2 = max(reshape(obj.U(1).x(2:Nx+1,:,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U(2).x(2:Nx+1,:,:),N,1), 0).*dMupxPh2;
                vecY1 = min(reshape(obj.U(1).y(:,1:Ny,:),N,1), 0).*dMupyPh1 + min(reshape(obj.U(2).y(:,1:Ny,:),N,1), 0).*dMupyPh2;
                vecY2 = max(reshape(obj.U(1).y(:,2:Ny+1,:),N,1), 0).*dMupyPh1 + max(reshape(obj.U(2).y(:,2:Ny+1,:),N,1), 0).*dMupyPh2;
                vecZ1 = min(reshape(obj.U(1).z(:,:,1:Nz),N,1), 0).*dMupzPh1 + min(reshape(obj.U(2).z(:,:,1:Nz),N,1), 0).*dMupzPh2;
                vecZ2 = max(reshape(obj.U(1).z(:,:,2:Nz+1),N,1), 0).*dMupzPh1 + max(reshape(obj.U(2).z(:,:,2:Nz+1),N,1), 0).*dMupzPh2;
                acc = pv/dt .* (obj.drhoTdp .* z(:,i));
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2+vecX2-vecY1-vecX1-vecZ1+acc, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                Jp{i} = Jp{i} + spdiags(DiagVecs, DiagIndx, N, N);
                
                %% 2. Component i composition block
                for j=1:obj.NofComponents-1
                    dMupxPh1 = obj.UpWind(1).x*...
                        ( obj.dMob(:, 1, j) .* x(:,(i-1)*2+1) .* rho(:,1)...
                        + obj.Mob(:, 1) .* obj.dxdz(:,(i-1)*2+1, j) .* rho(:,1)...
                        + obj.Mob(:, 1) .* x(:, (i-1)*2+1) .* obj.drhodz(:, 1));
                    dMupyPh1 = obj.UpWind(1).y*...
                        ( obj.dMob(:, 1, j) .* x(:, (i-1)*2+1) .* rho(:,1)...
                        + obj.Mob(:, 1) .* obj.dxdz(:, (i-1)*2+1, j) .* rho(:,1)...
                        + obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodz(:, 1));
                    dMupzPh1 = obj.UpWind(1).z*...
                        ( obj.dMob(:, 1, j) .* x(:,(i-1)*2+1) .* rho(:,1)...
                        + obj.Mob(:,1) .* obj.dxdz(:,(i-1)*2+1) .* rho(:,1)...
                        + obj.Mob(:,1) .* x(:,(i-1)*2+1) .* obj.drhodz(:, 1));
                    dMupxPh2 = obj.UpWind(2).x*...
                        ( obj.dMob(:, 2, j) .* x(:,(i-1)*2+2) .* rho(:,2)...
                        + obj.Mob(:,2) .* obj.dxdz(:,(i-1)*2+2) .* rho(:,2)...
                        + obj.Mob(:,2) .* x(:,(i-1)*2+2) .* obj.drhodz(:, 2));
                    dMupyPh2 = obj.UpWind(2).y*...
                        ( obj.dMob(:, 2, j) .* x(:,(i-1)*2+2) .* rho(:,2)...
                        + obj.Mob(:, 2) .* obj.dxdz(:,(i-1)*2+2, j) .* rho(:,2)...
                        + obj.Mob(:, 2) .* x(:, (i-1)*2+2) .* obj.drhodz(:, 2));
                    dMupzPh2 = obj.UpWind(2).z*...
                        ( obj.dMob(:, 2, j) .* x(:, (i-1)*2+2) .* rho(:,2)...
                        + obj.Mob(:,2) .* obj.dxdz(:, (i-1)*2+2, j) .* rho(:,2)...
                        + obj.Mob(:,2) .* x(:, (i-1)*2+2) .* obj.drhodz(:, 2));
                    
                    vecX1 = min(reshape(obj.U(1).x(1:Nx,:,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U(2).x(1:Nx,:,:),N,1), 0).*dMupxPh2;
                    vecX2 = max(reshape(obj.U(1).x(2:Nx+1,:,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U(2).x(2:Nx+1,:,:),N,1), 0).*dMupxPh2;
                    vecY1 = min(reshape(obj.U(1).y(:,1:Ny,:),N,1), 0).*dMupyPh1 + min(reshape(obj.U(2).y(:,1:Ny,:),N,1), 0).*dMupyPh2;
                    vecY2 = max(reshape(obj.U(1).y(:,2:Ny+1,:),N,1), 0).*dMupyPh1 + max(reshape(obj.U(2).y(:,2:Ny+1,:),N,1), 0).*dMupyPh2;
                    vecZ1 = min(reshape(obj.U(1).z(:,:,1:Nz),N,1), 0).*dMupzPh1 + min(reshape(obj.U(2).z(:,:,1:Nz),N,1), 0).*dMupzPh2;
                    vecZ2 = max(reshape(obj.U(1).z(:,:,2:Nz+1),N,1), 0).*dMupzPh1 + max(reshape(obj.U(2).z(:,:,2:Nz+1),N,1), 0).*dMupzPh2;
                    
                    if i == obj.NofComponents
                        acc = - pv/dt .* rhoT + pv/dt .* z(:,i) .* obj.drhoTdz(:, j);
                    elseif i == j
                        acc = pv/dt .* rhoT + pv/dt .* z(:,i) .* obj.drhoTdz(:, j);
                    end
                    
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    Jz{i,j} = spdiags(DiagVecs,DiagIndx, N, N);
                    % Capillarity
                    Jz{i,j} = Jz{i,j} - obj.Tph{i,1} * spdiags(obj.dPc, 0, N, N);
                end
            end
            
            %% 3. Add wells to each block
            [Jp, Jz] = obj.AddWellsToJacobian(Jp, Jz, P(:,obj.NofPhases), rho, x, ProductionSystem.Wells, ProductionSystem.Reservoir.K);
            
            %% Full Jacobian
            if obj.NofComponents <= 2 
            Jacobian = [Jp{1}, Jz{1,1}
                        Jp{2}, Jz{2,1}];
            else
            % for now let's consider a max of 3 components!! 
            Jacobian = [Jp{1}, Jz{1,1}, Jz{1,2};...
                        Jp{2}, Jz{2,1}, Jz{2,2};...
                        Jp{3}, Jz{3,1}, Jz{3,2}];
            end
        end
        function [Jp, Jz] = AddWellsToJacobian(obj, Jp, Jz, p, rho, x, Wells, K)
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:Wells.NofInj
                a = Inj(i).Cells;
                for j=1:length(a)
                    for c=1:obj.NofComponents
                    Jp{c}(a(j),a(j)) = Jp{c}(a(j),a(j)) ...
                                       + Inj(i).PI * K(a(j)) * (Inj(i).Mob(:,1) * Inj(i).rho(j,1) * Inj(i).x(:,(c-1)*2+1) ...
                                       + Inj(i).Mob(:,2) *  Inj(i).rho(j,2) * Inj(i).x(:,(c-1)*2+2));
                    end
                end
            end
            
            %Producers
            for i=1:Wells.NofProd
                b = Prod(i).Cells;
                for j=1:length(b)
                    for c=1:obj.NofComponents
                        %Pressure blocks
                        Jp{c}(b(j),b(j)) = Jp{c}(b(j),b(j))...
                            + Prod(i).PI * K(b(j)) *...
                            ( obj.Mob(b(j), 1) * rho(b(j),1) * x(b(j), (c-1)*2+1) ...
                            + obj.Mob(b(j), 2) * rho(b(j),2) * x(b(j), (c-1)*2+2)...
                            )...
                            - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * x(b(j), (c-1)*2+1) * obj.drhodp(b(j),1) * (Prod(i).p(j) - p(b(j))) ...
                            - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * obj.dxdp(b(j), (c-1)*2+1) * rho(b(j),1) * (Prod(i).p(j) - p(b(j))) ...
                            - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * x(b(j), (c-1)*2+2) * obj.drhodp(b(j),2) * (Prod(i).p(j) - p(b(j)))...
                            - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * obj.dxdp(b(j), (c-1)*2+2) * rho(b(j),2) * (Prod(i).p(j) - p(b(j)));
                        for zI = 1:obj.NofComponents-1
                            Jz{c, zI}(b(j),b(j)) = Jz{c, zI}(b(j),b(j))...
                                - Prod(i).PI * K(b(j)) *...
                                ( obj.dMob(b(j), 1, zI) * rho(b(j), 1) * x(b(j), (c-1)*2+1) ...
                                + obj.dMob(b(j), 2, zI) * rho(b(j), 2) * x(b(j), (c-1)*2+2)...
                                + obj.Mob(b(j), 1) * rho(b(j), 1) * obj.dxdz(b(j), (c-1)*2+1, zI) ...
                                + obj.Mob(b(j), 2) * rho(b(j), 2) * obj.dxdz(b(j), (c-1)*2+2, zI)...
                                + obj.Mob(b(j), 1, zI) * obj.drhodz(b(j), 1) * x(b(j), (c-1)*2+1) ...
                                + obj.Mob(b(j), 2, zI) * obj.drhodz(b(j), 2) * x(b(j), (c-1)*2+2)...
                                ) * (Prod(i).p(j) - p(b(j)));
                        end
                    end
                end          
            end
        end
        function UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                Nm =  DiscretizationModel.ReservoirGrid.N;
                %% 1. Update matrix
                % 1.a Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(obj.NofPhases)]);
                Pm.update(delta(1:Nm));
                % 1.b Update z
                DeltaLast = zeros(Nm, 1);
                for c = 1:obj.NofComponents-1
                    Zm = ProductionSystem.Reservoir.State.Properties(['z_', num2str(c)]);
                    Zm.update(delta(c*Nm + 1:(c+1)*Nm));
                    Zm.Value = max(Zm.Value, 0);
                    Zm.Value = min(Zm.Value, 1);
                    DeltaLast = DeltaLast + delta(c*Nm + 1:(c+1)*Nm);
                end
                Zm = ProductionSystem.Reservoir.State.Properties(['z_', num2str(obj.NofPhases)]);
                Zm.update(-DeltaLast);
                Zm.Value = max(Zm.Value, 0);
                Zm.Value = min(Zm.Value, 1);
                % Update Remaining properties (x, S, rho, rhoT, Pc)
                obj.UpdatePhaseCompositions(ProductionSystem.Reservoir.State, FluidModel);
                
                %% 2. Update fractures pressure and densities
                if ProductionSystem.FracturesNetwork.Active
                    for i=1:ProductionSystem.FracturesNetwork.NofFractures
                        % 2.a Update Pressure
                        Pf = ProductionSystem.FracturesNetwork.Fractures(i).State.Properties(['P_', num2str(obj.NofPhases)]);
                        Pf.update(delta);
                        % 2.b Update z
                        DeltaLast = zeros(Nf(i), 1);
                        for ph = 1:obj.NofPhases-1
                            Zf = ProductionSystem.FracturesNetwork.Fractures(i).State.Properties(['z_', num2str(ph)]);
                            Zf.update(delta(c*Nf(i) + 1:(c+1)*Nf(i)));
                            % Remove values that are not physical
                            Zf.Value = max(Zf.Value, 0);
                            Zf.Value = min(Zf.Value, 1);
                            DeltaLast = DeltaLast + delta(c*Nf(i) + 1:(c+1)*Nf(i));
                        end
                        Zf = ProductionSystem.FracturesNetwork.Fractures(i).State.Properties(['z_', num2str(obj.NofPhases)]);
                        Zf.update(-DeltaLast);
                        % Remove values that are not physical
                        Zf.Value = max(Zf.Value, 0);
                        Zf.Value = min(Zf.Value, 1);
                        % Update Remaining properties (x, S, rho, rhoT, Pc)
                        obj.UpdatePhaseCompositions(ProductionSystem.FracturesNetwork.Fractures(i).State, FluidModel);
                    end
                end
            end
        end
    end
end