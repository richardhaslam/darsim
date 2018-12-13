% LTS Linear solver class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_linear_solver < handle
    properties
        Name
        Tol
        Maxit
        Iter
        CondNumber
        %Preconditioner
    end
    methods
        function obj = LTS_linear_solver(name, tol, maxit)
            obj.Name = name;
            obj.Tol = tol;
            obj.Maxit = maxit;
        end
        function SetUp(obj, ProductionSystem, DiscretizationModel, Residual)
            %obj.Preconditioner.Setup(ProductionSystem, DiscretizationModel);
        end
        function x = Solve(obj, A, rhs, ActCells)
            
            A = A(ActCells, ActCells);
            rhs = rhs(ActCells);
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
                    [xred, flag, relres, obj.Iter] = gmres(A, rhs, [], obj.Tol, obj.Maxit, L, U);
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
                    [xred, flag, relres, obj.Iter] = bicg(A, rhs, obj.Tol, obj.Maxit, L, U);
                case('direct')
                    flag = 0;
                    xred = A\rhs;
                otherwise
                    error('unsupported linear solver type');
            end
            
            
            x = zeros(size(ActCells));
            x(ActCells) = xred;
            
            disp(['LS time:', num2str(toc(start))]);
            if flag ~= 0
                disp(['WARNING: Linear solver did not converge. The residual norm is ', num2str(relres)]);
            end
        end
    end
end