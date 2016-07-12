% NL solver base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef NL_Solver < handle
properties
    MaxChops
    MaxIter
    Tol
    LinearSolver
end
properties (Access = private)
    
end
methods
    function Solve(obj, ProductionSystem, Formulation, dt, Summary)
        % Initialise objects
        Converged=0;
        CompConverged = 0;
        chops=0;
        while ((Converged==0 || CompConverged == 0) && chops <= obj.MaxChops)
            
            %Update fluid prowperties
                     
            % Compute residual
            [Residual] = Formulation.BuildResidual(ProductionSystem);
            
            %Print some info to the screen
            if (chops > 0)
                disp('Maximum number of iterations was reached: time-step was chopped');
                disp(['Restart Newton loop dt = ', num2str(dt)]);
            end
            disp(['Initial residual norm: ', num2str(norm(Residual, inf))]);
            disp('');
            disp('        ||Residual||   ||delta p||   ||delta S||');
            
            %Linear Solver Setup
            obj.LinearSolver.setup(ProductionSystem, Discretization);
            
            % NEWTON LOOP
            % Initialise Timers
            TimerConstruct = zeros(obj.MaxIter,1);
            TimerSolve = zeros(obj.MaxIter, 1);
            TimerInner = zeros(obj.MaxIter, 1);
            itCount = 1;
            while ((Converged==0 || CompConverged == 0)  && (itCount <= obj.MaxIter))
                
                % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
                start1 = tic;
                J = Formulation.BuildJacobian();
                TimerConstruct(itCount) = toc(start1);
                
                % 2. Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
                start2 = tic;
                [Delta, Delta_c] = obj.LinearSolver.Solve(J, Residual, N, ADM);
                TimerSolve(itCount) = toc(start2);
                
                % 3. Update Solution
                ProductionSystem = Formulation.UpdateSolution(ProductionSystem);
                
                % 5. Compute residual
                [Residual] = Formulation.BuildResidual(ProductionSystem);
                
                % 6. Check convergence criteria
                Converged = NewtonConvergence(itCount, Residual, Delta+Delta2, Status.p, Tol, N, ADM, Delta_c);
                itCount = itCount+1;
            end
            if (Converged == 0 || CompConverged == 0)
                dt = dt/2;
                chops = chops + 1;
            end
        end
        
        %Choose next time-step size
        if itCount < 4
            dtnext = 2*dt;
        elseif itCount > 8
            dtnext = dt/2;
        else
            dtnext = dt;
        end
        
        %Compute Injection and production fluxes for Injection  and Production curves
        for i=1:length(Inj)
            c = Inj(i).cells;
            switch (Inj(i).type)
                case('RateConstrained')
                    Inj(i).qw = sum (Mw(c)./(Mw(c)+Mnw(c)).*Inj(i).q);
                    Inj(i).qnw = sum (Inj(i).q - Inj(i).qw(c));
                case('PressureConstrained')
                    %Phases
                    Inj(i).qw =   sum(Inj(i).Mw.* Inj(i).rho(c,1).* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24;
                    Inj(i).qnw =  sum(Inj(i).Mo(i).* Inj(i).rho(c,2) .* Inj(i).PI .* Kvector(c).* (Inj(i).p - Status.p(c)))*3600*24;
                    %Components
                    Inj(i).qz1 = sum(Inj(i).x1(1) .* Inj(i).Mw .* Inj.rho(1) .* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24 +...
                        sum(Inj(i).x1(2) .* Inj(i).Mo(c) .* Inj(i).rho(2) .* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24;
                    Inj(i).qz2 = sum((1 - Inj(i).x1(1)) .* Inj(i).Mw .* Inj(i).rho(1) .* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24 +...
                        sum((1 - Inj(i).x1(2)) .* Inj(i).Mo .* Inj(i).rho(2) .* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24;
            end
        end
        
        for i=1:length(Prod)
            c = Prod(i).cells;
            switch (Prod(i).type)
                case('RateConstrained')
                    Prod(i).qw = sum (Mw(c)./(Mw(c)+Mnw(c)).*Prod(i).q);
                    Prod(i).qnw = sum (Prod(i).q - Prod(i).qw(c));
                case('PressureConstrained')
                    %Phases
                    Prod(i).qw =   sum(Mw(c).* Status.rho(c,1).* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c)));
                    Prod(i).qnw =  sum(Mnw(c).* Status.rho(c,2) .* Prod(i).PI .* Kvector(c).* (Prod(i).p - Status.p(c)));
                    %Components
                    Prod(i).qz1 = sum(Status.x1(c, 1) .* Mw(c) .* Status.rho(c,1) .* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c))) +...
                        sum(Status.x1(c, 2) .* Mnw(c) .* Status.rho(c,2) .* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c)));
                    Prod(i).qz2 = sum((1 - Status.x1(c, 1)) .* Mw(c) .* Status.rho(c,1) .* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c))) +...
                        sum((1 - Status.x1(c, 2)) .* Mnw(c) .* Status.rho(c,2) .* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c)));
            end
        end
        
        %% Stats, timers and Injection/Production data
        Summary.CouplingStats.SaveStats(Ndt, itCount-1, chops);
        Summary.CouplingStats.SaveTimers(Ndt, TimerConstruct, TimerSolve, TimerInner);
        Summary.SaveWellsData(Ndt+1, Inj, Prod, dt);
    end
end
end