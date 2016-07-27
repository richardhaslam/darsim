% ADM Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver_ADM < linear_solver
    properties
        maxLevel
        Prolp
        Prols
        Rest
    end
    methods
        function obj = linear_solver_ADM()
            obj.Rest = struct('matrix', {});
            obj.Prolp = struct('matrix', {});
            obj.Prols = struct('matrix', {});
        end
        function SetUp(obj, ProductionSystem, DiscretizationModel)
            % Choose where to coarsen and build ADM grid
            DiscretizationModel.AdaptGrid();
            % Construct R & P based on ADM grid
            [obj.Rest, obj.Prolp, obj.Prols] = DiscretizationModel.ConstructOperators();
            obj.maxLevel = DiscretizationModel.ADMGrid.level(end);         
        end
        function xf = Solve(obj, A, rhs)
            rhs_c = obj.RestrictResidual(rhs);
            [A_c] = obj.RestrictSystem(A);
            obj.x = -A_c\rhs_c;
            xf = obj.Prolong();
        end
        function Residual_c = RestrictResidual(obj, Residual)
            Residual_cp = obj.Rest(1).matrix * Residual(1:length(Residual)/2);
            Residual_cs = obj.Rest(1).matrix * Residual(length(Residual)/2+1:length(Residual));
            for i=2:obj.maxLevel
                Residual_cp = obj.Rest(i).matrix*Residual_cp;
                Residual_cs = obj.Rest(i).matrix*Residual_cs;
            end
            Residual_c = [Residual_cp; Residual_cs];
        end
        function J_c = RestrictSystem(obj, J)
            N = length(J)/2;  % for now it should work
            % Restrict system to DLGR Grid
            Jop = J(1:N, 1:N);
            Jos = J(1:N, N+1:2*N);
            Jwp = J(N+1:2*N, 1:N);
            Jws = J(N+1:2*N, N+1:2*N);
            Jopc = R(1).matrix*Jop*Pp(1).matrix;
            Josc = R(1).matrix*Jos*Ps(1).matrix;
            Jwpc = R(1).matrix*Jwp*Pp(1).matrix;
            Jwsc = R(1).matrix*Jws*Ps(1).matrix;
            for i=2:obj.maxLevel
                Jopc = obj.Rest(i).matrix*Jopc*obj.Prolp(i).matrix;
                Josc = obj.Rest(i).matrix*Josc*obj.Prols(i).matrix;
                Jwpc = obj.Rest(i).matrix*Jwpc*obj.Prolp(i).matrix;
                Jwsc = obj.Rest(i).matrix*Jwsc*obj.Prols(i).matrix;
            end
            J_c = [Jopc, Josc; Jwpc, Jwsc];
        end
        function Delta = Prolong(obj, Delta_c)
            %Prolongs pressure Solution from DLGR grid to fine Grid
            Nc = length(Delta_c)/2;
            Deltap = Pp(obj.maxLevel).matrix * Delta_c(1:Nc);
            Deltas = Ps(obj.maxLevel).matrix * Delta_c(Nc+1:2*Nc);
            for i = obj.maxLevel-1:-1:1
                Deltap = obj.Prolp(i).matrix * Deltap;
                Deltas = obj.Prols(i).matrix * Deltas;
            end
            Delta = [Deltap; Deltas];
        end
    end
end