% ADM Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_linear_solver_ADM < linear_solver_ADM
    properties  
        Rred
        Pred
    end
    methods
        function obj = LTS_linear_solver_ADM(name, tol, maxit)
            obj@linear_solver_ADM(name, tol, maxit);
        end
        
        function SetUp(obj, ProductionSystem, DiscretizationModel, Residual)
            % Select grid and get ADM Operators 
            [obj.R, obj.P] = obj.OperatorsAssembler.Assemble(DiscretizationModel, ProductionSystem, Residual);
            
            obj.Rred = obj.R;
            obj.Pred = obj.P;
            % Modify permeability field (for upscaling option)
            if obj.DLGR
                DiscretizationModel.ModifyPerm(ProductionSystem); 
            end            
            
            % Display the number of active grids
            fprintf('Number of ADM Active Grids: %1.0f (%2.2f Percent of nodes)\n', ...
                    DiscretizationModel.ADMGrid.Ntot, 100*DiscretizationModel.ADMGrid.Ntot/DiscretizationModel.GlobalGrids(1).N);
        end
        
        function LTS_SetUp(obj, DiscretizationModel, ActCells, l)
            %modify the restiction and prolongation operators in accord to
            %the active componets
            maxLevel = DiscretizationModel.maxLevel;
            ActCells = logical(ActCells);
            
            % select the levels to compute inside R
            grid_l  =  maxLevel - l;
            obj.Rred = obj.R(DiscretizationModel.ADMGrid.level == 0, ActCells);
            for subL = 1:grid_l
                obj.Rred = cat(1, obj.Rred, obj.R(DiscretizationModel.ADMGrid.level == (subL), ActCells));
            end
            obj.Rred = obj.Rred(any( obj.Rred, 2),:);
            obj.Pred = obj.Rred';
        end
            
        function xf = Solve(obj, A, rhs, ActCells)
            ActCells = logical(ActCells);
            A = A(ActCells, ActCells);
            rhs = rhs(ActCells);
            % Restrict system
            rhs_c = obj.Rred * rhs;
            A_c = obj.Rred * A * obj.Pred;
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
            xp = obj.Pred * x;
            xf = zeros(size(ActCells));
            xf(ActCells) = xp;
            % xf = A\rhs;
            
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