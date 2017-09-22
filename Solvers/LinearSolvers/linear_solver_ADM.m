% ADM Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 21 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver_ADM < linear_solver
    properties  
        R
        P
        Smooth = 0
        OperatorsAssembler
    end
    methods
        function obj = linear_solver_ADM(name, tol, maxit)
            obj@linear_solver(name, tol, maxit);
        end
        function SetUp(obj, DiscretizationModel)
            % Get ADM Operators
            [obj.R, obj.P] = obj.OperatorsAssembler.Assemble(DiscretizationModel.OperatorsHandler.ADMRest, DiscretizationModel.OperatorsHandler.ADMProl); 
        end
        function xf = Solve(obj, A, rhs)
            % Restrict system
            rhs_c = obj.R * rhs;
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

            % Prolong to fine-scale resolution
            xf = obj.P * x;
            
            % Should be ilu(0) smoothing but it's not working
            if obj.Smooth
                setup.type = 'nofill';
                setup.milu = 'off';
                setup.droptol = 0.1;
                [L, U] = ilu(A, setup);
                for i=1:10
                    %rhs_s = rhs - A*xf;
                    [xf, flag, relres, n_iter] = bicgstab(A, rhs, 1e-1, [], L, U, xf);
                    disp(n_iter);
                    %delta = U\(L\rhs_s);
                    %delta = A\rhs_s;
                    %xf = xf + delta;
                end
            end
        end
    end
end