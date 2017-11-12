% Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 22 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver < handle
    properties
        Name
        Tol
        Maxit
        Iter
        CondNumber
        %Preconditioner
    end
    methods
        function obj = linear_solver(name, tol, maxit)
            obj.Name = name;
            obj.Tol = tol;
            obj.Maxit = maxit;
        end
        function SetUp(obj, ProductionSystem, DiscretizationModel)
            %obj.Preconditioner.Setup(ProductionSystem, DiscretizationModel);
        end
        function x = Solve(obj, A, rhs)
            %obj.CondNumber = condest(A);
            obj.CondNumber = 1;
            start = tic;       
            switch (obj.Name)
                case('gmres')
                    % Set-up ilu preconditioner
                    setup.type = 'nofill';
                    setup.milu = 'off';
                    setup.droptol = 0.1;
                    Diagonal = diag(A);
                    i=find(Diagonal==0);
                    if ~isempty(i)
                        disp(i);
                    end
                    [L, U] = ilu(A, setup);
                    [x, flag, relres, obj.Iter] = gmres(A, rhs, [], obj.Tol, obj.Maxit, L, U);
                case('bicg')
                    % Set-up ilu preconditioner
                    setup.type = 'nofill';
                    setup.milu = 'off';
                    setup.droptol = 0.1;
                    Diagonal = diag(A);
                    i=find(Diagonal==0);
                    if ~isempty(i)
                        disp(i);
                    end
                    [L, U] = ilu(A, setup);
                    [x, flag, relres, obj.Iter] = bicg(A, rhs, obj.Tol, obj.Maxit, L, U);
                case('direct')
                    flag = 0;
                    x = A\rhs;
                otherwise
                    error('unsupported linear solver type');
            end
            disp(['LS time:', num2str(toc(start))]);
            if flag ~= 0
                disp(['WARNING: Linear solver did not converge. The residual norm is ', num2str(relres)]);
            end
        end
    end
end