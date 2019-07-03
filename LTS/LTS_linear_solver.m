% LTS Linear solver class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_linear_solver < linear_solver
    properties
       ActiveComponents 
    end
    methods
        function obj = LTS_linear_solver(name, tol, maxit)
            obj@linear_solver(name, tol, maxit);
        end
        function SetUp(obj, ProductionSystem, DiscretizationModel, ActCells, Nph)
            % I should set up the number of active components
            N = DiscretizationModel.N;
            obj.ActiveComponents = false(N*Nph, 1);
            End = 0;
            for ph = 1:Nph
                Start = End + 1;
                End = Start - 1 + N;  
                obj.ActiveComponents(Start:End) = ActCells;
            end
        end
        function x = Solve(obj, A, rhs, varargin)
            
            if length(varargin) == 1 
                obj.ActiveComponents = varargin{1};
            end
             
            A = A(obj.ActiveComponents, obj.ActiveComponents);
            rhs = rhs(obj.ActiveComponents);
            % obj.CondNumber = condest(A);
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
                    error('DARSim2 error: unsupported linear solver type');
            end
            
            
            x = zeros(size(obj.ActiveComponents));
            x(obj.ActiveComponents) = xred;
            
            disp(['LS time:', num2str(toc(start))]);
            if flag ~= 0
                disp(['WARNING: Linear solver did not converge. The residual norm is ', num2str(relres)]);
            end
        end
    end
end