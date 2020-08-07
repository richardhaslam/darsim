% Geothermal MultiPhase formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_MultiPhase_formulation < formulation
    properties
        MatrixAssembler
        
        dMobdp
%         dSdp
        dTdp
%         d2Td2p
        
        drhodh
        dMobdh
        dSdh
%         d2rhodh2
%         d2Mobdh2
        dTdh
%         d2Td2h
        
%         Kf % fluid thermal conductivity
        Tk % transmisibility of rock conductivity
        Th % transmisibility of rho .* h
        
        Thph
        Ghph
    end
    methods
        function obj = Geothermal_MultiPhase_formulation(NofPhases)
            obj@formulation();
            obj.Tph = cell(NofPhases,1);
            obj.Gph = cell(NofPhases,1);
            obj.MatrixAssembler = matrix_assembler_geothermal();
            % This is the matrix assembler for geothermal (multiphase)
        end
        function x = GetPrimaryUnknowns(obj, ProductionSystem, DiscretizationModel)
             Nt = DiscretizationModel.N; % total grid
             Nm = DiscretizationModel.ReservoirGrid.N; % matrix grid
             x = zeros(2 * Nt, 1); % 2 unknowns: P, H
             Index.Start = 1;
             Index.End = Nm;
                     
             % Get matrix unknowns
             x(Index.Start     :Index.End       ) = ProductionSystem.Reservoir.State.Properties('P_1').Value;
             x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.Reservoir.State.Properties('hTfluid').Value;
             
             % Get fracture unknowns
             for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                 Index.Start = Index.End+1;
                 Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                 x(Index.Start     :Index.End       ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_1').Value;
                 x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('hTfluid').Value;
             end
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            % 1. Computing properties
            obj.ComputeProperties(ProductionSystem, FluidModel);
            
            % 2. Reservoir Properties and Derivatives
            obj.ComputeDerivatives(ProductionSystem, FluidModel);
        end
        function ComputeProperties(obj, ProductionSystem, FluidModel)
            %% 1. Geothermal Properties

            %% 2. Reservoir Properties 
            % Reservoir
            FluidModel.GetPhaseEnthalpies(ProductionSystem.Reservoir.State);
            FluidModel.GetPhaseIndex(ProductionSystem.Reservoir.State);
            
            FluidModel.GetTemperature(ProductionSystem.Reservoir.State);
            FluidModel.GetPhaseViscosities(ProductionSystem.Reservoir.State); 
            FluidModel.GetPhaseDensities(ProductionSystem.Reservoir.State); 
            FluidModel.ComputePhaseSaturations(ProductionSystem.Reservoir.State); 
            FluidModel.correctEnthalpy(ProductionSystem.Reservoir.State);
            
            % Compute enthalpy of rock
            FluidModel.ComputeRockEnthalpy(ProductionSystem.Reservoir);
            % Compute Effective Thermal Conductivity coefficient
            FluidModel.ComputeThermalConductivity(ProductionSystem.Reservoir);
            
            % for phase mobilities; including rel.perm model (linear). Computed based on single saturation value as S = Sa+Sb
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State);
                        
            % Update Pc [THIS GUY UPDATES P_1; Pfffff....]
            FluidModel.ComputePc(ProductionSystem.Reservoir.State);
            
            %% 3. Fractures Properties
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                FluidModel.GetPhaseEnthalpies(ProductionSystem.FracturesNetwork.Fractures(f).State);
                FluidModel.GetPhaseIndex(ProductionSystem.FracturesNetwork.Fractures(f).State);
                
                FluidModel.GetTemperature(ProductionSystem.FracturesNetwork.Fractures(f).State);
                FluidModel.GetPhaseViscosities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                FluidModel.GetPhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                FluidModel.ComputePhaseSaturations(ProductionSystem.FracturesNetwork.Fractures(f).State);
                FluidModel.correctEnthalpy(ProductionSystem.FracturesNetwork.Fractures(f).State);
                
                % Compute enthalpy of rock
                FluidModel.ComputeRockEnthalpy(ProductionSystem.FracturesNetwork.Fractures(f));
                % Compute Effective Thermal Conductivity coefficient
                FluidModel.ComputeThermalConductivity(ProductionSystem.FracturesNetwork.Fractures(f));
                
                % for phase mobilities; including rel.perm model (linear). Computed based on single saturation value as S = Sa+Sb
                obj.Mob = [obj.Mob; FluidModel.ComputePhaseMobilities(ProductionSystem.FracturesNetwork.Fractures(f).State)];
                
                % Update Pc [THIS GUY UPDATES P_1; Pfffff....]
                FluidModel.ComputePc(ProductionSystem.FracturesNetwork.Fractures(f).State);
            end
        end
        function ComputeDerivatives(obj, ProductionSystem, FluidModel)
            % Pressure and Enthalpy indices for tables are updated when
            % updating the properties; if properties are not updated, there
            % is no point updating the derivatives
            
            % THIS WORKS FOR FRACTURES, BUT RUINS CONVERGENCE WHEN NO FRACTURES ARE INVOLVED
%             FluidModel.GetPhaseEnthalpies(ProductionSystem.Reservoir.State);
%             FluidModel.GetPhaseIndex(ProductionSystem.Reservoir.State);

            % Derivatives with respect to Pressure
            obj.drhodp = FluidModel.ComputeDrhoDp(ProductionSystem.Reservoir.State); 
            obj.dMobdp = FluidModel.ComputeDMobDp(ProductionSystem.Reservoir.State);    % re-evaluate
%             obj.dSdp = FluidModel.ComputeDSDp(ProductionSystem.Reservoir.State);
            obj.dTdp = FluidModel.ComputeDTDp(ProductionSystem.Reservoir.State);
%             obj.d2Td2p = FluidModel.ComputeD2TD2p(ProductionSystem.Reservoir.State);
            
            % Derivatives with respect to Enthalpy
            obj.drhodh = FluidModel.ComputeDrhoDh(ProductionSystem.Reservoir.State); 
            obj.dMobdh = FluidModel.ComputeDMobDh(ProductionSystem.Reservoir.State);   % re-evaluate
            obj.dSdh = FluidModel.ComputeDSDh(ProductionSystem.Reservoir.State);
            obj.dTdh = FluidModel.ComputeDTDh(ProductionSystem.Reservoir.State);
%             obj.d2Td2h = FluidModel.ComputeD2TD2h(ProductionSystem.Reservoir.State);
            
            
            %% 3. Fractures Derivatives
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                % THIS WORKS FOR FRACTURES, BUT RUINS CONVERGENCE WHEN NO FRACTURES ARE INVOLVED
%                 FluidModel.GetPhaseEnthalpies(ProductionSystem.FracturesNetwork.Fractures(f).State);
%                 FluidModel.GetPhaseIndex(ProductionSystem.FracturesNetwork.Fractures(f).State);

                obj.drhodp = [obj.drhodp;FluidModel.ComputeDrhoDp(ProductionSystem.FracturesNetwork.Fractures(f).State)];
                obj.dMobdp = [obj.dMobdp;FluidModel.ComputeDMobDp(ProductionSystem.FracturesNetwork.Fractures(f).State)];    % re-evaluate
                obj.dTdp = [obj.dTdp;FluidModel.ComputeDTDp(ProductionSystem.FracturesNetwork.Fractures(f).State)];
                
                % Derivatives with respect to Enthalpy
                obj.drhodh = [obj.drhodh;FluidModel.ComputeDrhoDh(ProductionSystem.FracturesNetwork.Fractures(f).State)];
                obj.dMobdh = [obj.dMobdh;FluidModel.ComputeDMobDh(ProductionSystem.FracturesNetwork.Fractures(f).State)];   % re-evaluate
                obj.dSdh = [obj.dSdh;FluidModel.ComputeDSDh(ProductionSystem.FracturesNetwork.Fractures(f).State)];
                obj.dTdh = [obj.dTdh;FluidModel.ComputeDTDh(ProductionSystem.FracturesNetwork.Fractures(f).State)];   
            end
        end
        function Residual_MB = BuildMediumFlowResidual(obj, Medium, Grid, dt, State0, Index, qw, qf, f)
            Index = Index.Start:Index.End;
            
            % Initialize residual of Mass Balance
            Residual_MB = zeros(Grid.N,1);
            
            % Create local variables from time step n
            P_old    = State0.Properties('P_2').Value(Index);
            
            % Create local variables from iteration step nu
            P_new    = Medium.State.Properties('P_2').Value;
            
            % Pore Volume & Rock Volume
            Medium.ComputePorosity(P_old);
            pv_old = Medium.Por  .* Grid.Volume;   % Old pore Volume
            
            Medium.ComputePorosity(P_new);
            pv_new = Medium.Por  .* Grid.Volume;   % New pore volume
            
            % Phase dependent fluxes
            Accumulation_fluid = zeros(Grid.N,1); % Initialize accumulation term
            Convection_Flux = zeros(Grid.N,1);   % Initialize convection flux
            SourceTerm_Flux = zeros(Grid.N,1);   % Initialize source term flux

            for ph=1:obj.NofPhases % Loop over the phases
                rho_old = State0.Properties(['rho_', num2str(ph)]).Value(Index);
                S_old = State0.Properties(['S_', num2str(ph)]).Value(Index);

                P_new = Medium.State.Properties(['P_', num2str(ph)]).Value;   % Phase pressure from iteration step nu
                rho_new = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                S_new = Medium.State.Properties(['S_', num2str(ph)]).Value;
                
                % Mass Accumulation Term for Fluid
                Accumulation_fluid = Accumulation_fluid + ( pv_new .* rho_new .* S_new - pv_old .* rho_old .* S_old ) / dt;
                
                % Thermal Convection flux (summation of all phases)
                Convection_Flux = Convection_Flux + obj.Tph{ph, 1+f} * P_new - obj.Gph{ph, 1+f} * Grid.Depth;
                
                % SourceTerm flux
                SourceTerm_Flux = SourceTerm_Flux + qw(Index, ph) + qf(Index, ph);
            end
            
            % Fill residual with all flux terms: LHS - RHS
            Residual_MB = Residual_MB + Accumulation_fluid + Convection_Flux - SourceTerm_Flux;
        
            % end of function
        end
        function Residual_EB = BuildMediumHeatResidual(obj, Medium, Grid, dt, State0, Index, qhw, qhf, RTf, f)
            Index = Index.Start:Index.End;
            
            % Initialize residual of Energy Balance
            Residual_EB = zeros(Grid.N,1);
            
            % Create local variables from time step n
            hRock_old   = State0.Properties('hRock').Value(Index);
            P_old       = State0.Properties('P_2').Value(Index);
            
            % Create local variables from iteration step nu
            hRock_new   = Medium.State.Properties('hRock').Value;
            P_new       = Medium.State.Properties('P_2').Value;

            % Pore Volume & Rock Volume
            Medium.ComputePorosity(P_old);
            mv_old = (1 - Medium.Por) .* Grid.Volume;   % Old rock volume
            pv_old =      Medium.Por  .* Grid.Volume;   % Old pore Volume
            
            Medium.ComputePorosity(P_new);
            mv_new = (1 - Medium.Por) .* Grid.Volume;   % New rock volume
            pv_new =      Medium.Por  .* Grid.Volume;   % New pore volume

            % Energy Accumulation Term for Rock
            Rho_rock = Medium.Rho;
            Accumulation_rock = ( mv_new .* Rho_rock .* hRock_new - mv_old .* Rho_rock .* hRock_old ) / dt;

            % Phase dependent fluxes
            Accumulation_fluid = zeros(Grid.N,1); % Initialize accumulation term
            Convection_Flux = zeros(Grid.N,1);    % Initialize convection flux
            Conduction_Flux = zeros(Grid.N,1);    % Initialize conduction flux
            SourceTerm_Flux = zeros(Grid.N,1);    % Initialize source term flux

            for ph=1:obj.NofPhases % Loop over the phases
                h_old   = State0.Properties(['h_',num2str(ph)]).Value(Index);
                rho_old = State0.Properties(['rho_',num2str(ph)]).Value(Index);
                S_old   = State0.Properties(['S_',num2str(ph)]).Value(Index);
                
                P_new = Medium.State.Properties(['P_', num2str(ph)]).Value;   % Phase pressure from iteration step nu
                hTfluid_new = Medium.State.Properties('hTfluid').Value;
                h_new = Medium.State.Properties(['h_',num2str(ph)]).Value;
                rho_new = Medium.State.Properties(['rho_',num2str(ph)]).Value;
                S_new = Medium.State.Properties(['S_',num2str(ph)]).Value;

                % Energy Accumulation Term for Flkuid
                Accumulation_fluid = Accumulation_fluid + ...
                    ( pv_new .* rho_new .* S_new .* h_new - pv_old .* rho_old .* S_old .* h_old ) / dt;
                
                % Thermal Convection flux (summation of all phases)
                Convection_Flux = Convection_Flux + obj.Thph{ph, 1+f} * P_new - obj.Ghph{ph, 1+f} * Grid.Depth;
                
                % SourceTerm flux
                SourceTerm_Flux = SourceTerm_Flux + qhw(Index, ph) + qhf(Index, ph);
            end
            
            % Thermal Conduction flux: CondEff = (1-phi)*D_rock + phi*S_water*D_water + phi*S_steam*D_steam 
            Conduction_Flux = Conduction_Flux + obj.Tk{1, 1+f} * ( obj.dTdp(Index) .* P_new + obj.dTdh(Index) .* hTfluid_new );
            
            % Fill residual with all flux terms: LHS - RHS
            Residual_EB = Residual_EB + Accumulation_rock + Accumulation_fluid + Conduction_Flux + Convection_Flux - SourceTerm_Flux; 
        
            % end of function
        end
        function ResidualFull = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [Qw, Qhw] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);              % Mass Flow flux between each two media
            Qhf= zeros(DiscretizationModel.N, obj.NofPhases);              % Heat Convection flux betweem each two media
            RTf= zeros(DiscretizationModel.N, 1);                          % Heat Conduction flux betweem each two media
            if ProductionSystem.FracturesNetwork.Active
                [Qf, Qhf, RTf] = obj.ComputeSourceTerms_frac_mat(ProductionSystem, DiscretizationModel);
            end
            RTf= zeros(DiscretizationModel.N, 1);

            % Initialise residual vector (2 * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
            else
                Nf = 0;
            end
            Nt = DiscretizationModel.N;
            ResidualFull = zeros( 2*Nt , 1 ); % We have only 2 equations now
            
            % Compute conductive tranmissibility
            obj.Tk{1,1} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.ReservoirGrid , ProductionSystem.Reservoir.State );
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                obj.Tk{1, 1+f} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.FracturesGrid.Grids(f) , ProductionSystem.FracturesNetwork.Fractures(f).State );
            end
            
            for ph=1:obj.NofPhases
                %% Transmissibilities of flow and heat
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                % Compute flow transmissibility and gravity
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
                % Compute (convective) heat transmissibility and gravity
                [obj.Thph{ph, 1}, obj.Ghph{ph, 1}] = obj.MatrixAssembler.ConvectiveHeatTransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + Nf(f) - 1;
                    % Compute flow transmissibility and gravity
                    [obj.Tph{ph, 1+f}, obj.Gph{ph, 1+f}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                        DiscretizationModel.FracturesGrid.Grids(f), ...
                        obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph), ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, ...
                        obj.GravityModel.RhoInt{ph, 1+f});
                    % Compute (convective) heat transmissibility and gravity
                    [obj.Thph{ph, 1+f}, obj.Ghph{ph, 1+f}] = obj.MatrixAssembler.ConvectiveHeatTransmissibilityMatrix( ...
                        DiscretizationModel.FracturesGrid.Grids(f), ...
                        obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph), ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['h_',num2str(ph)]).Value, ...
                        obj.GravityModel.RhoInt{ph, 1+f});
                end
            end
            %% Flow Residual
            % Reservoir
            Index.Start = 1;
            Index.End = Nm;
            Residualm = BuildMediumFlowResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qw, Qf, 0);
            ResidualFull(Index.Start: Index.End) = Residualm;
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End+1;
                Index.End = Index.Start + Nf(f) - 1;
                Residualf = BuildMediumFlowResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, Qw, Qf, f);
                ResidualFull(Index.Start:Index.End) = Residualf;
            end

            %% Heat Residual
            % Reservoir
            Index.Start = Index.End + 1;
            Index.End = Index.Start + Nm - 1;
            Index_r.Start = 1;
            Index_r.End = Nm;
            Residualm = BuildMediumHeatResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index_r, Qhw, Qhf, RTf, 0);
            ResidualFull(Index.Start: Index.End) = Residualm;
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End + 1;
                Index.End = Index.Start + Nf(f) - 1;
                Index_r.Start = Index_r.End + 1;
                Index_r.End = Index_r.Start + Nf(f) - 1;
                Residualf = BuildMediumHeatResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index_r, Qhw, Qhf, RTf, f);
                ResidualFull(Index.Start:Index.End) = Residualf;
            end
        end
        function [J_MB_P , J_MB_H] = BuildMediumFlowJacobian(obj, Medium, Wells, Grid, dt, Index, f)
            % Create local variables
            Index = Index.Start:Index.End;

            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            % Initialize Jacobian blocks for Mass Balance
            J_MB_P = sparse(N,N);
            J_MB_H = sparse(N,N);
            
            % Porosity and its derivative 
            P = Medium.State.Properties('P_2').Value;
            Medium.ComputeDerPorosity(P);
            phi = Medium.Por;
            dphidp = Medium.DPor;
                 
            % Phase dependent derivatives
            % Initialize derivatives
            Derivative_Accumulation_Fluid = zeros(N,1);
            Transmissibility = zeros(N,1);
            Derivative_Transmissibility = zeros(N,1);

            % Phase dependent derivatives
            for ph=1:obj.NofPhases % loop over all phases
                rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                S = Medium.State.Properties(['S_',num2str(ph)]).Value;
                
                %% J_MB_P Block
                % Derivative of Fluid Accumulation
                vec = (Grid.Volume/dt) .* ( ...
                                            dphidp .* rho                  .* S + ...
                                            phi    .* obj.drhodp(Index,ph) .* S ...
                                           );
                
                Derivative_Accumulation_Fluid = Derivative_Accumulation_Fluid + spdiags(vec, 0, N, N);

                % Derivative of Fluid Mass Convection Flux:
                Transmissibility = Transmissibility + obj.Tph{ph,1+f}; % (Transmissibility Term) x (Derivative of P_nu)

                % (Derivative of Transmissibility Term) x (P_nu)
                dMupx = obj.UpWind{ph,1+f}.x *( obj.Mob(Index, ph)    .* obj.drhodp(Index, ph) + ...
                                                obj.dMobdp(Index, ph) .* rho            );
                dMupy = obj.UpWind{ph,1+f}.y *( obj.Mob(Index, ph)    .* obj.drhodp(Index, ph) + ...
                                                obj.dMobdp(Index, ph) .* rho            );
                dMupz = obj.UpWind{ph,1+f}.z *( obj.Mob(Index, ph)    .* obj.drhodp(Index, ph) + ...
                                                obj.dMobdp(Index, ph) .* rho            );
                
                % Because of multiplication with velocity obj.U, we are multiplying it with grad(P^nu).
                vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;

                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                
                % Construct (Derivative of Transmissibility Term)
                Derivative_Transmissibility = Derivative_Transmissibility + spdiags(DiagVecs, DiagIndx, N, N);
                
            end
            
            % Construction of J_MB_P
            J_MB_P = J_MB_P + Derivative_Accumulation_Fluid + Transmissibility + Derivative_Transmissibility;
                
            
            %% J_MB_H Block
            
            % Phase dependent derivatives
            % Initialize derivatives
            Derivative_Accumulation_Fluid = zeros(N,1);
            Transmissibility = zeros(N,1);
            Derivative_Transmissibility = zeros(N,1);
    
            for ph=1:obj.NofPhases % loop over all phases
                rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                S = Medium.State.Properties(['S_',num2str(ph)]).Value;

                % Derivative of Fluid Accumulation
                vec = (Grid.Volume/dt) .* ( ...
                                            phi .* obj.drhodh(Index,ph) .* S + ...
                                            phi .* rho                  .* obj.dSdh(Index,ph) ...
                                           );
                Derivative_Accumulation_Fluid = Derivative_Accumulation_Fluid + spdiags(vec,0,N,N);
                % *** we could use rhoT here, and take the Derivative_Accumulation_Fluid out of the phase loop ***
                
                % Derivative of Fluid Mass Convection Flux:
                % (Transmissibility Term) x (Derivative of P_nu)
                Transmissibility = Transmissibility + 0;

                % (Derivative of Transmissibility Term) x (P_nu)
                dMupx = obj.UpWind{ph,1+f}.x * ( obj.Mob(Index, ph)    .* obj.drhodh(Index, ph) + ...
                                                 obj.dMobdh(Index, ph) .* rho            );
                dMupy = obj.UpWind{ph,1+f}.y * ( obj.Mob(Index, ph)    .* obj.drhodh(Index, ph) + ...
                                                 obj.dMobdh(Index, ph) .* rho            );
                dMupz = obj.UpWind{ph,1+f}.z * ( obj.Mob(Index, ph)    .* obj.drhodh(Index, ph) + ...
                                                 obj.dMobdh(Index, ph) .* rho            );
                
                % Because of multiplication with velocity obj.U, we are multiplying it with grad(P^nu).
                vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
                
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                
                % Construct (Derivative of Transmissibility Term)
                Derivative_Transmissibility = Derivative_Transmissibility + spdiags(DiagVecs, DiagIndx, N, N);
                
            end % end phase loop
            
            % construction of J_MB_H
            J_MB_H = J_MB_H + Derivative_Accumulation_Fluid + Transmissibility + Derivative_Transmissibility;
            
                            
            % Add Source Term to all Jacobian blocks
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [J_MB_P, J_MB_H] = obj.AddWellsToMassBalanceJacobian(J_MB_P, J_MB_H, Medium.State, Wells, Medium.K(:,1));
            end
            
            % end of function
        end
        function [J_EB_P , J_EB_H] = BuildMediumHeatJacobian(obj, Medium, Wells, Grid, dt, Index, f)
            % Create local variables
            Index = Index.Start:Index.End;

            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            % Initialize Jacobian blocks for Energy Balance
            J_EB_P = sparse(N,N);
            J_EB_H = sparse(N,N);
            
            % Porosity and its derivative 
            phi = Medium.Por;
            Medium.ComputeDerPorosity(Medium.State.Properties('P_2').Value);
            dphidp = Medium.DPor;
            
            %% J_EB_P Block
            % Derivative of Rock Energy Accumulation
            hRock = Medium.State.Properties('hRock').Value; 
            rhoRock = Medium.Rho;
            
            vec = (Grid.Volume/dt) .* (-1) .* dphidp .* rhoRock .* hRock; 
            Derivative_Accumulation_Rock = spdiags(vec, 0, N, N);

            % Phase dependent derivatives
            % Initialize derivatives
            Derivative_Accumulation_Fluid = zeros(N,1);
            Convective_Transmissibility = zeros(N,1);
            Derivative_Convective_Transmissibility = zeros(N,1);
            Conductive_Transmissibility = zeros(N,1);
            Derivative_Conductive_Transmissibility = zeros(N,1);
            
            % Total Fluid Enthalpy 
            hTfluid = Medium.State.Properties('hTfluid').Value;
            P = Medium.State.Properties('P_2').Value;
            
            for ph=1:obj.NofPhases % loop over all phases
                rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                h = Medium.State.Properties(['h_',num2str(ph)]).Value;
                S = Medium.State.Properties(['S_',num2str(ph)]).Value;

                % Derivative of Fluid Energy Accumulation
                vec = (Grid.Volume/dt) .* ( ...
                                            dphidp .* rho              .* S              .* h + ... % should be phase enthalpy h
                                            phi    .* obj.drhodp(Index,ph) .* S              .* h   ...
                                           );

%                                             phi    .* rho              .* obj.dSdp(Index,ph) .* h ...

                
                Derivative_Accumulation_Fluid = Derivative_Accumulation_Fluid + spdiags(vec,0,N,N);

                % Derivative of Fluid Energy Convection Flux:
                % (Transmissibility Term) x (Derivative of P_nu)
                Convective_Transmissibility = Convective_Transmissibility + obj.Thph{ph, 1+f};

                % (Derivative of Transmissibility Term) x (P_nu)
                dMupx = obj.UpWind{ph,1+f}.x*( obj.Mob(Index, ph)    .* obj.drhodp(Index, ph) .* h + ...
                                               obj.dMobdp(Index, ph) .* rho                   .* h  );   % ADD DhDp AS h IS FUNCTION OF Pressure...??                                               
                dMupy = obj.UpWind{ph,1+f}.y*( obj.Mob(Index, ph)    .* obj.drhodp(Index, ph) .* h + ...
                                               obj.dMobdp(Index, ph) .* rho                   .* h  );
                dMupz = obj.UpWind{ph,1+f}.z*( obj.Mob(Index, ph)    .* obj.drhodp(Index, ph) .* h + ...
                                               obj.dMobdp(Index, ph) .* rho                   .* h  );
                
                % Because of multiplication with velocity obj.U, we are multiplying it with grad(P^nu).
                vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;

                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                
                % Construct (Derivative of Transmissibility Term)
                Derivative_Convective_Transmissibility = Derivative_Convective_Transmissibility + spdiags(DiagVecs, DiagIndx, N, N);
                                
            end % end of phase loop
               
            % Derivative of Fluid Energy Convection Flux:
            % (Transmissibility Term) x (Derivative of P_nu)
            Conductive_Transmissibility = Conductive_Transmissibility + obj.Tk{1,1+f} .* obj.dTdp(Index); 
            
            % (Derivative of Transmissibility Term) x (P_nu)
            Derivative_Conductive_Transmissibility = Derivative_Conductive_Transmissibility ;%+ obj.Tk{1,1+f} .* obj.d2Td2p .* P;
            % should this P be delta(P)? Same as velocity?

            % construction of J_EB_P 
            J_EB_P = J_EB_P + Derivative_Accumulation_Rock + Derivative_Accumulation_Fluid + Convective_Transmissibility + ...
                              Derivative_Convective_Transmissibility + Conductive_Transmissibility + Derivative_Conductive_Transmissibility;

                      
            %% J_EB_H Block
            % Derivative of Rock Energy Accumulation
            vec = zeros(N,1); 
            Derivative_Accumulation_Rock = spdiags(vec, 0, N, N);

            % Phase dependent derivatives
            % Initialize derivatives
            Derivative_Accumulation_Fluid = zeros(N,1);
            Convective_Transmissibility = zeros(N,1);
            Derivative_Convective_Transmissibility = zeros(N,1);
            Conductive_Transmissibility = zeros(N,1);
            Derivative_Conductive_Transmissibility = zeros(N,1);
            
            % Total Fluid Enthalpy 
            hTfluid = Medium.State.Properties('hTfluid').Value;
            
            for ph=1:obj.NofPhases % loop over all phases
                rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                h = Medium.State.Properties(['h_',num2str(ph)]).Value;
                S = Medium.State.Properties(['S_',num2str(ph)]).Value;

                % Derivative of Fluid Energy Accumulation
                vec = (Grid.Volume/dt) .* ( ...
                                            phi .* obj.drhodh(Index,ph) .* S              .* h + ...
                                            phi .* rho              .* obj.dSdh(Index,ph) .* h + ...
                                            phi .* rho              .* S              .* 1         ... 
                                           );           
                % *** This is the part where phi*rho*S was added to the accumulation term TWICE ***
                
                Derivative_Accumulation_Fluid = Derivative_Accumulation_Fluid + spdiags(vec,0,N,N);

                % Derivative of Fluid Energy Convection Flux:
                % (Transmissibility Term) x (Derivative of P_nu)
                Convective_Transmissibility = Convective_Transmissibility + 0;

                % (Derivative of Transmissibility Term) x (P_nu)
                dMupx = obj.UpWind{ph,1+f}.x*( obj.Mob(Index, ph)    .* obj.drhodh(Index, ph) .* h + ...
                                               obj.dMobdh(Index, ph) .* rho                   .* h + ...
                                               obj.Mob(Index, ph)    .* rho                   .* 1 );
                dMupy = obj.UpWind{ph,1+f}.y*( obj.Mob(Index, ph)    .* obj.drhodh(Index, ph) .* h + ...
                                               obj.dMobdh(Index, ph) .* rho                   .* h + ...
                                               obj.Mob(Index, ph)    .* rho                   .* 1 );
                dMupz = obj.UpWind{ph,1+f}.z*( obj.Mob(Index, ph)    .* obj.drhodh(Index, ph) .* h + ...
                                               obj.dMobdh(Index, ph) .* rho                   .* h + ...
                                               obj.Mob(Index, ph)    .* rho                   .* 1 );

                % Because of multiplication with velocity obj.U, we are multiplying it with grad(P^nu).
                vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;

                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                
                % Construct (Derivative of Transmissibility Term)
                Derivative_Convective_Transmissibility = Derivative_Convective_Transmissibility + spdiags(DiagVecs, DiagIndx, N, N);
                                
            end % end of phase loop
               
            % Derivative of Fluid Energy Convection Flux:
            % (Transmissibility Term) x (Derivative of P_nu)
            Conductive_Transmissibility = Conductive_Transmissibility + obj.Tk{1,1+f} .* obj.dTdh(Index); 
            
            % (Derivative of Transmissibility Term) x (P_nu)
            Derivative_Conductive_Transmissibility = Derivative_Conductive_Transmissibility ;%+ obj.Tk{1,1+f} .* obj.d2Td2h .* hTfluid;

            % construction of J_EB_P 
            J_EB_H = J_EB_H + Derivative_Accumulation_Rock + Derivative_Accumulation_Fluid + Convective_Transmissibility + ...
                              Derivative_Convective_Transmissibility + Conductive_Transmissibility + Derivative_Conductive_Transmissibility;
                                        
            % Add Source Term to all Jacobian blocks
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [J_EB_P, J_EB_H] = obj.AddWellsToEnergyBalanceJacobian(J_EB_P, J_EB_H, Medium.State, Wells, Medium.K(:,1));
            end
            
            % end of function
        end
        function JacobianFull = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %% Jacobian's assembly for MultiPhase geothermal reservoir
            % | JPP  JPT | dP |    | RP |
            % | JTP  JTT | dT | = -| RT |
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            Reservoir = ProductionSystem.Reservoir;
            Fractures = ProductionSystem.FracturesNetwork.Fractures;
            Wells = ProductionSystem.Wells;
            
            %% Jacobian of the reservoir
            Index.Start = 1;
            Index.End = Nm;
            % Flow Jacobian blocks
            [J_MB_P, J_MB_H] = BuildMediumFlowJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir
            % Heat Jacobian blocks
            [J_EB_P, J_EB_H] = BuildMediumHeatJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir
            
            %% Jacobian of the fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Nf = DiscretizationModel.FracturesGrid.N;
                Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                Index.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, Nf(f));
                % Flow Jacobian Block
                [Jf_MB_P, Jf_MB_H] = BuildMediumFlowJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f);
                J_MB_P  = blkdiag(J_MB_P, Jf_MB_P);
                J_MB_H  = blkdiag(J_MB_H, Jf_MB_H);
                % Heat Jacobian Block
                [Jf_EB_P, Jf_EB_H] = BuildMediumHeatJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f);
                J_EB_P  = blkdiag(J_EB_P, Jf_EB_P);
                J_EB_H  = blkdiag(J_EB_H, Jf_EB_H);
            end
            
            %% ADD frac-matrix and frac-frac connections in Jacobian blocks
            % Global variables
            if ProductionSystem.FracturesNetwork.Active
                FineGrid = [DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids];
            else
                FineGrid = DiscretizationModel.ReservoirGrid;
            end
            P = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'P_'); % useful for cross connections assembly
            rho = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'rho_'); % useful for cross connections assembly
            h = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'h_'); % useful for cross connections assembly
            for c = 1:length(DiscretizationModel.CrossConnections)
                T_Geo = DiscretizationModel.CrossConnections(c).T_Geo; % geological transmissibility; = harmonic average of permeability
                T_Geo_Cond = DiscretizationModel.CrossConnections(c).T_Geo_Cond;
                UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                i = c + Nm;
                j = DiscretizationModel.CrossConnections(c).Cells;
                
                for ph = 1:obj.NofPhases
                    %% J_MB_P Coupling
                    % qf = T_Geo * Upwind * Mob * rho * (dP)
                    
                    % fracture - matrix coupling; qf in matrix equation
                    J_MB_P_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( rho(j,ph) + obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) ) );
                    % matrix - fracture coupling; qf in fracture equation
                    J_MB_P_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( rho(i,ph) + obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_MB_P_conn  = J_MB_P_1_conn + J_MB_P_2_conn;
                    
                    % frac - mat or frac1 - frac2
%                     if sum(J_MB_P(i, j))~=0 || sum(J_MB_P(j, i))~=0
%                         error('J_PP(i, j) or J_PP(j, i) is not zero!');
%                     end
                    J_MB_P(i, j) = J_MB_P(i, j) + J_MB_P_conn';
                    J_MB_P(i, i) = J_MB_P(i, i) - sum(J_MB_P_conn);
                    % mat-frac or frac2 - frac1
                    J_MB_P(j, i) = J_MB_P(j, i) + J_MB_P_conn;
                    J_MB_P(sub2ind([Nt, Nt], j, j)) = J_MB_P(sub2ind([Nt, Nt], j, j)) - J_MB_P_conn;
                    
                    %% J_MB_H Coupling
                    % qf = T_Geo * Upwind * Mob * rho * (dP)
                    
                    J_MB_H_1_conn = - T_Geo .*  UpWind(:, ph) .* (P(j, ph) - P(i, ph)) .* ( rho(j, ph)        .* obj.dMobdh(j, ph) + ...
                                                                                            obj.drhodh(j, ph) .* obj.Mob(j, ph)     );
                    J_MB_H_2_conn = - T_Geo .* ~UpWind(:, ph) .* (P(j, ph) - P(i, ph)) .* ( rho(i, ph)        .* obj.dMobdh(i, ph) + ...
                                                                                            obj.drhodh(i, ph) .* obj.Mob(i, ph)     );
                    
                    % frac - mat or frac1 - frac2
%                     if sum(J_MB_H(i, j))~=0 || sum(J_MB_H(j, i))~=0
%                         error('J_PT(i, j) or J_PT(j, i) is not zero!');
%                     end
                    J_MB_H(i, j) = J_MB_H(i, j) + J_MB_H_1_conn';
                    J_MB_H(i, i) = J_MB_H(i, i) + sum(J_MB_H_2_conn);
                    % mat-frac or frac2 - frac1
                    J_MB_H(j, i) = J_MB_H(j, i) - J_MB_H_2_conn;
                    % diag of mat or frac2
                    J_MB_H(sub2ind([Nt, Nt], j, j)) = J_MB_H(sub2ind([Nt, Nt], j, j)) - J_MB_H_1_conn;
                    
                    %% J_EB_P Coupling
                    % qf = T_Geo * Upwind * Mob * rho * h * (dp)
                    
                    J_EB_P_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* h(j, ph) .* ( rho(j,ph) .* 1                            + ...
                                                                                                 obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_EB_P_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* h(i, ph) .* ( rho(i,ph) .* 1                            + ...
                                                                                                 obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_EB_P_conn = J_EB_P_1_conn + J_EB_P_2_conn;
                    % frac - mat or frac1 - frac2
%                     if sum(J_EB_P(i, j))~=0 || sum(J_EB_P(j, i))~=0
%                         error('J_TP(i, j) or J_TP(j, i) is not zero!');
%                     end
                    J_EB_P(i, j) = J_EB_P(i, j) + J_EB_P_conn';
                    J_EB_P(i, i) = J_EB_P(i, i) - sum(J_EB_P_conn);
                    % mat-frac or frac2 - frac1
                    J_EB_P(j, i) = J_EB_P(j, i) + J_EB_P_conn;
                    % diag of mat or frac2
                    J_EB_P(sub2ind([Nt, Nt], j, j)) = J_EB_P(sub2ind([Nt, Nt], j, j)) - J_EB_P_conn;
                    
                    %% J_EB_H Coupling
                    % qf = T_Geo * Upwind * Mob * rho * h * (dp)
                    
                    J_EB_H_1_conn = - T_Geo .*  UpWind(:, ph) .* (P(j, ph) - P(i, ph)) .* ( obj.dMobdh(j, ph) .* rho(j, ph)        .* h(j, ph) + ...
                                                                                            obj.Mob(j, ph)    .* obj.drhodh(j, ph) .* h(j, ph) + ...
                                                                                            obj.Mob(j, ph)    .* rho(j, ph)        .* 1         );
                    J_EB_H_2_conn = - T_Geo .* ~UpWind(:, ph) .* (P(j, ph) - P(i, ph)) .* ( obj.dMobdh(i, ph) .* rho(i, ph)        .* h(i, ph) + ...
                                                                                            obj.Mob(i, ph)    .* obj.drhodh(i, ph) .* h(i, ph) + ...
                                                                                            obj.Mob(i, ph)    .* rho(i, ph)        .* 1         );
                    % frac - mat or frac1 - frac2
%                     if sum(J_EB_H(i, j))~=0 || sum(J_EB_H(j, i))~=0
%                         error('J_TT(i, j) or J_TT(j, i) is not zero!');
%                     end
                    J_EB_H(i, j) = J_EB_H(i, j) + J_EB_H_1_conn';
                    J_EB_H(i, i) = J_EB_H(i, i) + sum(J_EB_H_2_conn);
                    % mat-frac or frac2 - frac1
                    J_EB_H(j, i) = J_EB_H(j, i) - J_EB_H_2_conn;
                    % diag of mat or frac2
                    J_EB_H(sub2ind([Nt, Nt], j, j)) = J_EB_H(sub2ind([Nt, Nt], j, j)) - J_EB_H_1_conn;
                end
            end
            % Build & Stack Jacobian
            JacobianFull = [J_MB_P, J_MB_H ; J_EB_P, J_EB_H];
        end
        function ConstrainedPressureResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [Qw, ~] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);
            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            for ph=1:obj.NofPhases
                %% Transmissibilities of flow and heat
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                % Compute flow transmissibility and gravity
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
            end
                
            %% Mass Balance Residual
            Index.Start = 1;
            Index.End = Nm;
            R_MB = BuildMediumFlowResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qw, Qf, 0);
            
            %% Mass Balance Jacobian (J_MB_P)
            Nm = DiscretizationModel.ReservoirGrid.N;
            Reservoir = ProductionSystem.Reservoir;
            Wells = ProductionSystem.Wells;
            Index.Start = 1;
            Index.End = Nm;
            % Flow Jacobian blocks
            [J_MB_P, ~] = BuildMediumFlowJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir

            % Solve for deltaP
            deltaP = J_MB_P\(-R_MB);
                        
            %% Update Reservoir State
            Pm = ProductionSystem.Reservoir.State.Properties('P_2');
            Pm.update(deltaP(1:Nm));
            % Update ALL properties
            obj.ComputeProperties(ProductionSystem, FluidModel)            
            % 7. Update wells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            
%             fprintf('Norm CPR = %1.5e \n', norm(R_MB) );
%             fprintf('Norm Delta P = %1.5e \n', norm(deltaP./max(Pm.Value)) );

        end
        function ConstrainedEnthalpyResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [~,Qhw] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qhf = zeros(DiscretizationModel.N, obj.NofPhases);
            RTf = zeros(DiscretizationModel.N, obj.NofPhases);
            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            % Compute conductive tranmissibility
%             obj.Tk{1,1} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.ReservoirGrid );
            obj.Tk{1,1} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.ReservoirGrid , ProductionSystem.Reservoir.State );
            
            for ph=1:obj.NofPhases
                %% Transmissibilities of flow and heat
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                % Compute (convective) heat transmissibility and gravity
                [obj.Thph{ph, 1}, obj.Ghph{ph, 1}] = obj.MatrixAssembler.ConvectiveHeatTransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
            end
            %% Energy Balance Residual
            % Reservoir
            Index.Start = 1;
            Index.End = Nm;
            R_EB = BuildMediumHeatResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qhw, Qhf, RTf, 0);
                            
            %% Mass Balance Jacobian (J_MB_P)
            Nm = DiscretizationModel.ReservoirGrid.N;
            Reservoir = ProductionSystem.Reservoir;
            Wells = ProductionSystem.Wells;
            Index.Start = 1;
            Index.End = Nm;
            % Heat Jacobian blocks
            [~, J_EB_H] = BuildMediumHeatJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir

            % Solve for deltaH
            deltaH = J_EB_H\(-R_EB);
            
            %% Update Reservoir State
            Hm = ProductionSystem.Reservoir.State.Properties('hTfluid');
            Hm.update(deltaH(1:Nm));
            % Update ALL properties
            obj.ComputeProperties(ProductionSystem, FluidModel)            
            % 7. Update wells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end       
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                % Inflexion point correction
%                 delta = obj.UpdateEnthalpy(delta, ProductionSystem, FluidModel, DiscretizationModel);

                Nm = DiscretizationModel.ReservoirGrid.N;
                Nt = DiscretizationModel.N;
                deltaP = delta(1:Nt);
                deltaH = delta(Nt+1:2*Nt);
                
                %% Update reservoir state
                % 1. Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties('P_2');
                Pm.update(deltaP(1:Nm));
                % 2. Update Enthalpy
                Hm = ProductionSystem.Reservoir.State.Properties('hTfluid');
                Hm.update(deltaH(1:Nm));
                                
                %% Update fracture state
                % 1. Update Pressure
                if ProductionSystem.FracturesNetwork.Active
                    EP = Nm;
                    Nf = DiscretizationModel.FracturesGrid.N;
                    for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                        IP = EP+1;
                        EP = IP + Nf(f) - 1;
                        Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_2');
                        Pf.update(deltaP(IP:EP));
                        % 2. Update Fluid Temperature
                        H = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('hTfluid');
                        H.update(deltaH(IP:EP));                  
                    end
                end
                
                % Update properties and derivatives
                obj.ComputeProperties(ProductionSystem, FluidModel)            
                % 7. Update wells
                ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            end
        end
        function [Qw, Qhw]= ComputeSourceTerms(obj, N, Wells)
            Qw = zeros(N, obj.NofPhases);
            Qhw = zeros(N, obj.NofPhases);  
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj(i).Cells;
                Qw(c, :) = Wells.Inj(i).QPhases(:,:); % MB source term
                Qhw(c, :) = Wells.Inj(i).Qh(:,:);     % EB source term 
                % function is UpdateState in
                % ProductionSystem>Wells>injector_pressure.m [class]
            end
            %Producers
            for i=1:Wells.NofProd
                c = Wells.Prod(i).Cells;
                Qw(c, :) = Wells.Prod(i).QPhases(:,:);
                Qhw(c, :) = Wells.Prod(i).Qh(:,:);
            end
        end
        function [Qf, Qhf, RTf]= ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel)
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);              % Mass Flow flux between each two media 
            Qhf = zeros(DiscretizationModel.N, obj.NofPhases);             % Heat Convection flux betweem each two media
            RTf = zeros(DiscretizationModel.N, 1);                         % Heat Conduction flux betweem each two media
            Nm = DiscretizationModel.ReservoirGrid.N;

            % Global variables
            P = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'P_'); % useful for cross connections assembly
            rho = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'rho_'); % useful for cross connections assembly
            h = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'h_'); % useful for cross connections assembly
            T = ProductionSystem.CreateGlobalSinglePhaseVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], 'T'); % useful for cross connections assembly
            for c=1:length(DiscretizationModel.CrossConnections)
                j = DiscretizationModel.CrossConnections(c).Cells;
                i = c + Nm;
                T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                T_Geo_Cond = DiscretizationModel.CrossConnections(c).T_Geo_Cond;
                UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                for ph=1:obj.NofPhases
                    Qf(i, ph) = Qf(i, ph) + sum(T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* (P(j, ph) - P(i, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* (P(j, ph) - P(i, ph))) );
                    Qf(j, ph) = Qf(j, ph) + T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* (P(i, ph) - P(j, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* (P(i, ph) - P(j, ph)) );
                    
                    Qhf(i, ph) = Qhf(i, ph) + sum(T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* h(j,ph) .* (P(j, ph) - P(i, ph)) + ...
                        ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* h(i,ph) .* (P(j, ph) - P(i, ph))) );
                    Qhf(j, ph) = Qhf(j, ph) + T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* h(j,ph) .* (P(i, ph) - P(j, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* h(i,ph) .* (P(i, ph) - P(j, ph)) );
                end
                RTf(i) = RTf(i) + sum( T_Geo_Cond .* (T(j) - T(i)) );
                RTf(j) = RTf(j) + T_Geo_Cond .* (T(i) - T(j));
            end
        end
        function [J_MB_P, J_MB_H] = AddWellsToMassBalanceJacobian(obj, J_MB_P, J_MB_H, State, Wells, K)
            % Define Local handles
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                dQdp = Inj(i).ComputeWellMassFluxDerivativeWithRespectToPressure(K, obj.NofPhases);
                dQdh = Inj(i).ComputeWellMassFluxDerivativeWithRespectToEnthalpy(obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(a)
                        % add derivative of inj well to J_MB_P
                        J_MB_P(a(j),a(j)) = J_MB_P(a(j),a(j)) - dQdp(j, ph);
                        J_MB_H(a(j),a(j)) = J_MB_H(a(j),a(j)) - dQdh(j, ph);
                    end
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                dQdp = Prod(i).ComputeWellMassFluxDerivativeWithRespectToPressure(State, K, obj.Mob, obj.drhodp, obj.dMobdp, obj.NofPhases);
                dQdh = Prod(i).ComputeWellMassFluxDerivativeWithRespectToEnthalpy(State, K, obj.Mob, obj.drhodh, obj.dMobdh, obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(b)
                        % add derivative of prod well to J_MB_P & J_MB_H
                        J_MB_P(b(j),b(j)) = J_MB_P(b(j),b(j)) - dQdp(j, ph);
                        J_MB_H(b(j),b(j)) = J_MB_H(b(j),b(j)) - dQdh(j, ph);
                    end
                end
            end
        end
        function [J_EB_P, J_EB_H] = AddWellsToEnergyBalanceJacobian(obj, J_EB_P, J_EB_H, State, Wells, K)
            % Define Local handles
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                dQhdp = Inj(i).ComputeWellHeatFluxDerivativeWithRespectToPressure(State, K, obj.NofPhases);
                dQhdh = Inj(i).ComputeWellHeatFluxDerivativeWithRespectToEnthalpy(State, K, obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(a)
                        % add derivative of inj well to J_EB_P
                        J_EB_P(a(j),a(j)) = J_EB_P(a(j),a(j)) - dQhdp(j, ph);
                        J_EB_H(a(j),a(j)) = J_EB_H(a(j),a(j)) - dQhdh(j, ph);
                    end
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                dQhdp = Prod(i).ComputeWellHeatFluxDerivativeWithRespectToPressure(State, K, obj.Mob, obj.drhodp, obj.dMobdp, obj.NofPhases);
                dQhdh = Prod(i).ComputeWellHeatFluxDerivativeWithRespectToEnthalpy(State, K, obj.Mob, obj.drhodh, obj.dMobdh, obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(b)
                        % add derivative of prod well to J_EB_P & J_EB_H
                        J_EB_P(b(j),b(j)) = J_EB_P(b(j),b(j)) - dQhdp(j, ph);
                        J_EB_H(b(j),b(j)) = J_EB_H(b(j),b(j)) - dQhdh(j, ph);
                    end
                end
            end
        end
        function d2Fdh2 = ComputeSecondDerivativeHeatConvectionFlux(obj, Medium, FluidModel, hTfluid)
            d2Fdh2 = zeros(size(obj.Mob,1),1);
            
            % Construct P and H grid for interp functions 
            [Hgrid,Pgrid] = meshgrid(FluidModel.Htable,FluidModel.Ptable);


            STable = FluidModel.TablePH.('S_1'); 
            S1 = FluidModel.Phases(1).GetSaturation(Pgrid, Hgrid, STable, hTfluid, Medium.State.Properties('P_2').Value); 
            kr = FluidModel.RelPermModel.ComputeRelPerm(FluidModel.Phases, S1);

            for i=1:obj.NofPhases
                % Density
                rhoTable = FluidModel.TablePH.(['rho_', num2str(i)]);
                rho(:,i) = FluidModel.Phases(i).GetDensity(Pgrid, Hgrid, rhoTable, hTfluid, Medium.State.Properties('P_2').Value);
                drhodh(:,i) = FluidModel.Phases(i).ComputeDrhoDh(Pgrid, Hgrid, rhoTable, hTfluid, Medium.State.Properties('P_2').Value);
                
                % Phase Enthalpy 
                PhaseEnthalpyTable = FluidModel.TablePH.(['H_',num2str(i)]);
                h(:,i) = FluidModel.Phases(i).GetPhaseEnthalpy(FluidModel.Ptable, PhaseEnthalpyTable, Medium.State.Properties('P_2').Value);    

                % 2nd derivative density wrt enthalpy
                d2rhodh2(:,i) = FluidModel.Phases(i).ComputeD2rhoDh2(Pgrid, Hgrid, rhoTable, hTfluid, Medium.State.Properties('P_2').Value);
                
                % Viscosity
                muTable = FluidModel.TablePH.(['mu_',num2str(i)]);
                mu(:,i) = FluidModel.Phases(i).GetViscosity(Pgrid, Hgrid, muTable, hTfluid, Medium.State.Properties('P_2').Value);
                
                % Mobility
                Mob(:,i) = kr(:,i)./mu(:,i);

                % 1st derivatives Viscosity and Mobility
                dmudh(:,i) = FluidModel.Phases(i).ComputeDmuDh(Pgrid, Hgrid, muTable, hTfluid, Medium.State.Properties('P_2').Value);
                dMobdh(:,i) = ( -1 .* dmudh(:,i) .* kr(:,i) ) ./ mu(:,i).^2;
                
                % 2nd derivatives Viscosity and Mobility
                d2mudh2(:,i) = FluidModel.Phases(i).ComputeD2muDh2(Pgrid, Hgrid, muTable, hTfluid, Medium.State.Properties('P_2').Value);
                d2Mobdh2(:,i) = -1 .* kr(:,i) .* ( d2mudh2(:,i) ./ mu(:,i).^2 - 2 .* dmudh(:,i) ./ mu(:,i).^3 );   
            end
            
            % Correction to Mobility and its derivatives
            Mob(isnan(Mob))=0; 
            dMobdh(isnan(dMobdh))=0;
            d2Mobdh2(isnan(d2Mobdh2))=0;

            % F = rho .* h .* Mob;
            % dFdh = drhodh .* h .* Mob + rho .* h .* dMobdh;
            
            % d2Fd2h = d2rhod2h .* h .* Mob + drhodh .* h .* dMobdh + drhodh .* h .* dMobdh + rho .* h .* d2Mobd2h;
            % per phase
            for ph=1:obj.NofPhases
                d2Fdh2 = d2Fdh2 + ...
                         d2rhodh2(:,ph) .* h(:,ph) .* Mob(:,ph) + ...
                         drhodh(:,ph)   .* h(:,ph) .* dMobdh(:,ph) + ...
                         drhodh(:,ph)   .* h(:,ph) .* dMobdh(:,ph) + ...
                         rho(:,ph)      .* h(:,ph) .* d2Mobdh2(:,ph) ;            
            end 
        end
        function delta = UpdateEnthalpy(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            % Detection of inflexion point (you have a bunch of zeros...)
            Nt = DiscretizationModel.N;
            deltaP_Old = delta(1:Nt);
            deltaH_Old = delta(Nt+1:2*Nt);

            H_old = ProductionSystem.Reservoir.State.Properties('hTfluid').Value;
            H_new = H_old + deltaH_Old;
            
            % Flux Correction - Mousa and Me
            d2Fdh2_old = obj.ComputeSecondDerivativeHeatConvectionFlux( ProductionSystem.Reservoir, FluidModel, H_old);
            d2Fdh2_new = obj.ComputeSecondDerivativeHeatConvectionFlux( ProductionSystem.Reservoir, FluidModel, H_new);
            H_new = H_new.*(d2Fdh2_new.*d2Fdh2_old >= 0) + 0.5*(H_new + H_old).*(d2Fdh2_new.*d2Fdh2_old < 0);
            deltaH_New = H_new - H_old;
            
            
            % This is the actual update of the delta solution (guess)
            deltaP_New = deltaP_Old .* (deltaH_New./deltaH_Old); 
            % Not sure about this part (line 823); it appears that the pressure is/must also be corrected,
            % i.e. steered into the same direction as enthalpy (formerly temperature)
            delta(   1:  Nt) = deltaP_New;
            delta(Nt+1:2*Nt) = deltaH_New;
            
        end
        function ComputeTotalFluxes(obj, ProductionSystem, DiscretizationModel)
            % this is virtual call
        end
        function AverageMassOnCoarseBlocks(obj, ProductionSystem, FineGrid, FluidModel, ADMRest)
            T = ProductionSystem.CreateGlobalSinglePhaseVariables(FineGrid, 'T'); % No phase index for this variable
            rho = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'rho_');
            % Perform Average for ADM
            T_rest = ADMRest * T;
            rho_rest = ADMRest * rho;     
            T_rho_rest = ADMRest * (T .* rho);  
            Tav = ADMRest' * ( T_rest ./ sum(ADMRest, 2));
            delta_T = Tav - T;
                     
            rho_av = ADMRest' * ( rho_rest ./ sum(ADMRest, 2));
            delta_rho = rho_av - rho;

            % Updating the variables in reservoir
            Start=1;
            End = FineGrid(1).N;
            T = ProductionSystem.Reservoir.State.Properties('T');
            rho = ProductionSystem.Reservoir.State.Properties('rho_1');
            T.update(delta_T(Start:End));
            rho.update(delta_rho(Start:End));
            
            % Updating the variables in fractures
            for frac = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                Start = End + 1;
                End = Start + FineGrid(frac+1).N - 1;
                T = ProductionSystem.FracturesNetwork.Fractures(frac).State.Properties('T');
                rho = ProductionSystem.FracturesNetwork.Fractures(frac).State.Properties('rho_1');
                T.update(delta_T(Start:End));
                rho.update(delta_rho(Start:End));
            end
        end
        function CFL = ComputeCFLNumber(obj, ProductionSystem, DiscretizationModel, dt)
            CFL = 0;
        end
    end
end