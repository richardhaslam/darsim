% Iterative Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 November 2016
%Last modified: 16 November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver_iterative < linear_solver
    properties
        Name
        Tol
        Maxit
    end
    methods
        function obj = linear_solver_iterative(name, tol, maxit)
            obj.Name = name;
            obj.Tol = tol;
            obj.Maxit = maxit;
        end
        function SetUp(obj, ProductionSystem, DiscretizationModel)
        end
        function x = Solve(obj, A, rhs)
            % Solve system with iterative linear solver
            %obj.CondNumber = condest(A);
            obj.CondNumber = 1;
            
            % Set-up preconditioner
            setup.type = 'nofill';
            setup.milu = 'off';
            setup.droptol = 0.1; 
            [L, U] = ilu(A, setup);
            switch (obj.Name)
                case('gmres')
                    restart = 3;
                    [x, flag, relres] = gmres(A, rhs, restart, obj.Tol, obj.Maxit, L, U);
                case('bicg')
                    [x, flag, relres] = bicg(A, rhs, obj.Tol, obj.Maxit, L, U);
            end
            if flag ~= 0
                disp(['Linear solver did not converge. The residual norm is ', num2str(relres)]);
            end
        end
    end
end