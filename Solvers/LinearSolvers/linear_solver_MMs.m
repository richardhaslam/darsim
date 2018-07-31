% MMs Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver_MMs < linear_solver
    properties
        R % Restriction operator
        P % Prolongation operator
        C % Correction functions operator
        OperatorsAssembler
        MSFE = 0;
    end
    methods
        function obj = linear_solver_MMs(name, tol, maxit)
            obj@linear_solver(name, tol, maxit);
        end
        function SetUp(obj, ProductionSystem, DiscretizationModel, Residual)
            if obj.MSFE
                obj.R = DiscretizationModel.OperatorsHandler.P';
            else
                obj.R = DiscretizationModel.OperatorsHandler.R;
            end
            obj.P = DiscretizationModel.OperatorsHandler.P;
            obj.C = DiscretizationModel.OperatorsHandler.C;
        end
        function xf = Solve(obj, A, rhs)
            % Restrict system
            % not using correction functions atm
            %if size(A,2) == size(obj.C,1)
             %   rhs_c = obj.R * (rhs - A * obj.C * rhs);
            %else
                rhs_c = obj.R * rhs;
            %end
            A_c = obj.R * A * obj.P;
 
            % Solve Coarse System
            %start = tic;
            switch (obj.Name)
                case('gmres')
                    % Set-up preconditioner
                    setup.type = 'nofill';
                    setup.milu = 'off';
                    setup.droptol = 0.1; 
                    [L, U] = ilu(A_c, setup);
                    [x, flag, relres, obj.Iter] = gmres(A_c, rhs_c, [], obj.Tol, min(obj.Maxit, size(A_c,1)), L, U);
                case('bicg')
                    % Set-up preconditioner
                    setup.type = 'nofill';
                    setup.milu = 'off';
                    setup.droptol = 0.1;
                    [L, U] = ilu(A_c, setup);
                    [x, flag, relres, obj.Iter] = bicg(A_c, rhs_c, obj.Tol, min(obj.Maxit, size(A_c,1)), L, U);
                otherwise
                    flag = 0;
                    x = A_c\rhs_c;
            end
            %disp(['Solving coarse system took ', num2str(toc(start)), ' s']);
            if flag == 1
                disp(['WARNING: Linear solver did not converge. The residual norm is ', num2str(relres)]);
            end
            
            % Prolong to fine-scale resolution and apply correction
            % functions
            %if size(obj.C,2) == size(rhs,1)
                %xf = obj.P * x + obj.C * rhs;
            %else
                xf = obj.P * x;
            %end
        end
    end
end