%  Natural variable Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef NaturalVar_formulation < fim_formulation
    properties
        T1
        T2
        Tw
        Tnw
        Mob
        dMob
        dPc
    end
    methods
        function ComputeDerivatives(obj, ProductionSystem, FluidModel)
            obj.dMob = FluidModel.MobilityDerivative(ProductionSystem.Reservoir.State.S);
        end
        function Residual = BuildResidual(obj, ProductionSystem, FluidModel, DiscretizationModel, dt, State0)
            
            %Initialise local variables
            s_old = State0.s;
            x1_old = State0.x1;
            x2_old = 1 - x1_old;
            p = ProductionSystem.Reservoir.State.p;
            s = ProductionSystem.Reservoir.State.S;
            x1 = ProductionSystem.Reservoir.State.x1;
            x2 = 1 - x1;
            s2 = 1 - s;
            s2_old = 1 - s_old;
            Pc = Status.pc;
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            %Density
            rho = State.rho;
            rho_old = State0.rho;
            
            %Accumulation term
            A = speye(Grid.N)*pv/dt;
            %Component 1
            m1 = x1(:,1) .* rho(:,1) .* s + x1(:,2) .* rho(:,2) .* s2;
            m1_old = x1_old(:,1) .* rho_old(:,1) .* s_old + x1_old(:,2) .* rho_old(:,2) .* s2_old;
            %Component 2
            m2 = x2(:,1) .* rho(:,1) .* s + x2(:,2) .* rho(:,2) .* s2;
            m2_old = x2_old(:,1) .* rho_old(:,1).* s_old + x2_old(:,2) .* rho_old(:,2) .* s2_old;
            
            %Convective term
            T1 = obj.TransmissibilityMatrix (Grid, obj.UpWindW, obj.UpWindNw, obj.Mob, rho, x1);
            T2 = obj.TransmissibilityMatrix (Grid, obj.UpWindW, obj.UpWindNw, obj.Mob, rho, x2);
            
            %Capillarity
            Tw = obj.TransmissibilityMatrix (Grid, obj.UpWindW, obj.UpWindNw, obj.Mob, rho, [ones(Grid.N,1), zeros(Grid.N,1)]);
            
            %Gravity
            G = ComputeGravityTerm(Grid.N);
            
            %Source terms
            [q1, q2] = ComputeWells(Grid.N, Inj, Prod, K, Status, Mnw, Mw, Rho);
            
            %% RESIDUAL
            %Component 1
            R1 = A * m1 - A * m1_old...           %Accumulation term
                + T1 * p...                       %Convective term
                - x1(:,1) .* (Tw*Pc)...           %Capillarity
                + G*p...                          %Gravity
                - q1;                             %Source terms
            %Component 2
            R2 = A * m2 - A * m2_old...
                + T2 * p...                       %Convective term
                - x2(:,1) .* (Tw*Pc)...           %Capillarity
                + G*p...                          %Gravity
                - q2;                             %Source terms
            
            %Stick them together
            Residual = [R1; R2];
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
        end
        function Status = UpdateState(obj, delta, Status, FluidModel)
            
        end
        function  T = TransmissibilityMatrix(obj, Grid, M, rho)
            %%%Transmissibility matrix construction
            Tx = zeros(Grid.Nx+1, Grid.Ny);
            Ty = zeros(Grid.Nx, Grid.Ny+1);
            
            %Apply upwind operator
            Mupx = UpWindW.x*(Mw .* Rho(:,1) .* x(:,1)) + UpWindNw.x*(Mnw .* Rho(:,2) .* x(:,2));
            Mupy = UpWindW.y*(Mw .* Rho(:,1) .* x(:,1)) + UpWindNw.y*(Mnw .* Rho(:,2) .* x(:,2));
            Mupx = reshape(Mupx, Grid.Nx, Grid.Ny);
            Mupy = reshape(Mupy, Grid.Nx, Grid.Ny);
            Tx(2:Grid.Nx,:)= Grid.Tx(2:Grid.Nx,:).*Mupx(1:Grid.Nx-1,:);
            Ty(:,2:Grid.Ny)= Grid.Ty(:,2:Grid.Ny).*Mupy(:,1:Grid.Ny-1);
            
            %Construct matrix
            x1 = reshape(Tx(1:Grid.Nx,:), Grid.N, 1);
            x2 = reshape(Tx(2:Grid.Nx+1,:), Grid.N, 1);
            y1 = reshape(Ty(:,1:Grid.Ny), Grid.N, 1);
            y2 = reshape(Ty(:,2:Grid.Ny+1), Grid.N, 1);
            DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
            DiagIndx = [-Grid.Nx, -1, 0, 1, Grid.Nx];
            T = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
        end
        function [q1, q2] = ComputeSourceTerms(obj, N, Wells)
        end
        function [Jwp, JwS, Jnwp, JnwS] = AddWellsToJacobian(obj, Jwp, JwS, Jnwp, JnwS, State, Wells, K)
        end
        function G = ComputeGravityTerm(obj, N)
        end
end