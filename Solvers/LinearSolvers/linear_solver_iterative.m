% Iterative Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 November 2016
%Last modified: 17 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver_iterative < linear_solver
    properties
        Name
        Tol
        Maxit
        Iter
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
            start = tic;
            Diagonal = diag(A);
            i=find(Diagonal==0);
            if ~isempty(i)
                disp(i);
            end
            [L, U] = ilu(A, setup);
            switch (obj.Name)
                case('gmres')
                    [x, flag, relres, obj.Iter] = gmres(A, rhs, [], obj.Tol, obj.Maxit, L, U);
                case('bicg')
                    [x, flag, relres, obj.Iter] = bicg(A, rhs, obj.Tol, obj.Maxit, L, U);
            end
            disp(['LS time:', num2str(toc(start))]);
            if flag ~= 0
                disp(['WARNING: Linear solver did not converge. The residual norm is ', num2str(relres)]);
            end
        end
    end
end