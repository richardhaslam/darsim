% Immiscible Fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 18 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_fluid_model < fluid_model
    properties
        name = 'Black Oil';
        Pref % Pref for Rs computation
        Pdim % Dimensionless pressure
        FlashSettings
    end
    methods
        function obj = BO_fluid_model(n_phases, n_comp)
            obj@fluid_model(n_phases, n_comp);
            obj.Pref = 1e7;
        end
        function AddFlash(obj, flash)
            obj.FlashSettings = flash;
        end
        function Status = InitializeReservoir(obj, Status)
            % Define initial values
            P_init = ones(length(Status.p), 1) * 5e7;
            z_init = ones(length(Status.p), 1)*0.1;
            
            % 1. Assign initial valus
            Status.p = Status.p .* P_init;
            Status.z(:,1) = z_init;
            Status.z(:,2) = 1 - z_init;
            
            
            % 2. Update Composition of the phases (Flash)
            SinglePhase = obj.Flash(Status);
            
            % 2.a Compute Phase Density
            for i=1:obj.NofPhases
                Status.rho(:, i) = obj.Phases(i).ComputeDensity(Status.p);
            end
            
            %3. Update S based on components mass balance
            obj.ComputePhaseSaturation(Status, SinglePhase);
            
            % 4. Total Density
            Status.rhoT = obj.ComputeTotalDensity(Status.S, Status.rho);
            
            % 5. Compute initial Pc
            Status.pc = obj.ComputePc(Status.S);
        end
        function InitializeInjectors(obj, Inj)
            % Loop over all injectors
            for i=1:length(Inj)
                Inj(i).z = [1 0];
                SinglePhase = obj.Flash(Inj(i));
                for ph=1:obj.NofPhases
                    Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i).p);
                end
                obj.ComputePhaseSaturation(Inj(i), SinglePhase);
                Inj(i).x2 = 1 - Inj(i).x1;
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).S);   
            end            
        end
        function SinglePhase = Flash(obj, Status)
            % Define SinglePhase objects
            SinglePhase.onlyliquid = zeros(length(Status.p), 1);
            SinglePhase.onlyvapor = zeros(length(Status.p), 1);
            
            %% - Dimensionless pressure!
            obj.Pdim = Status.p/obj.Pref;
            
            %% - Solve for x's
            Status.x1(:,2) = 1 - (800./(800 + 100*(0.2*obj.Pdim(:,1) + 0.2)));  %More pressure less comp 1 (oil) in phase 1 (oil) (actually more gas pushed in really)
            Status.x1(:,1) = 1 - 0;                                         %light component is all in gas
            
            %Recognize single phase cells and fix their xs to be equal to z
            SinglePhase.onlyliquid(Status.x1(:, 2) >= Status.z(:,1)) = 1;
            Status.x1(SinglePhase.onlyliquid == 1, 2) = Status.z(SinglePhase.onlyliquid == 1, 1);
        end
    end
end