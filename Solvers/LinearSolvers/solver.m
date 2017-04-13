classdef solver < handle
    properties
        name
        Tol
        Maxit
    end
    methods
        function x = Solve(obj, A, rhs)
            switch(obj.name)
                case('direct')
                    x = A\rhs;
                case('gmres')
                    [x, flag, relres, obj.Iter] = gmres(A, rhs, [], obj.Tol, obj.Maxit, L, U);
                case('bicg')
                    [x, flag, relres, obj.Iter] = bicg(A, rhs, obj.Tol, obj.Maxit, L, U);
            end
            if flag ~= 0
                disp(['WARNING: Linear solver did not converge. The residual norm is ', num2str(relres)]);
            end
        end
    end
end