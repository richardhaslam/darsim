%Add wells to Jacobean matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 13 June 2016
%Last Modified: 13 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J1p, J2p, J1S, J2S] = AddWellsToJacobeanComp(J1p, J2p, J1S, J2S, Inj, Prod, K, Status, Rho, dRho, Mw, Mnw, dMw, dMnw)
x1 = Status.x1;
x2 = 1 - Status.x1;
p = Status.p;

%Injectors
for i=1:length(Inj)
    a = Inj(i).cells;
    if strcmp(Inj(i).type, 'PressureConstrained')
        for j=1:length(a)
            J1p(a(j),a(j)) = J1p(a(j),a(j)) ...
                + Inj(i).PI * K(a(j)) * (Inj(i).Mw * Inj(i).Rho(1,1) * Inj(i).x1(1,1) ...
                + Inj(i).Mo *  Inj(i).Rho(1,2) * Inj(i).x1(1,2));
            J2p(a(j),a(j)) = J2p(a(j),a(j)) ...
                + Inj(i).PI * K(a(j)) * (Inj(i).Mw * Inj(i).Rho(1,1) * Inj(i).x2(1,1) ...
                + Inj(i).Mo *  Inj(i).Rho(1,2) * Inj(i).x2(1,2));
        end
    end
end

%Producers
for i=1:length(Prod)
    b = Prod(i).cells;
    switch(Prod(i).type)
        case('PressureConstrained')
            for j=1:length(b)
                %Pressure blocks
                J1p(b(j),b(j)) = J1p(b(j),b(j))...
                                 + Prod(i).PI * K(b(j)) * (Mw(b(j)) * Rho(b(j),1) * x1(b(j),1) ...
                                 + Mnw(b(j)) * Rho(b(j),2) * x1(b(j),2))...
                                 - Prod(i).PI * K(b(j)) * Mw(b(j)) * x1(b(j),1) * dRho(b(j),1) * (Prod(i).p - p(b(j))) ...
                                 - Prod(i).PI * K(b(j)) * Mnw(b(j)) * x1(b(j),2) * dRho(b(j),2) * (Prod(i).p - p(b(j)));
                J2p(b(j),b(j)) = J2p(b(j),b(j)) ...
                                 + Prod(i).PI * K(b(j)) * (Mw(b(j)) * Rho(b(j),1) * x2(b(j),1) ...
                                 + Mnw(b(j)) * Rho(b(j),2) * x2(b(j),2))...
                                 - Prod(i).PI * K(b(j)) * Mw(b(j)) * x2(b(j),1) * dRho(b(j),1) * (Prod(i).p - p(b(j))) ...
                                 - Prod(i).PI * K(b(j)) * Mnw(b(j)) * x2(b(j),2) * dRho(b(j),2) * (Prod(i).p - p(b(j)));
                %Saturation blocks
                J1S(b(j),b(j)) = J1S(b(j),b(j))...
                                 - Prod(i).PI * K(b(j)) * (dMw(b(j)) * Rho(b(j),1) * x1(b(j),1) ...
                                 + dMnw(b(j)) * Rho(b(j),2) * x1(b(j),2)) * (Prod(i).p - p(b(j)));
                J2S(b(j),b(j)) = J2S(b(j),b(j)) ...
                                 - Prod(i).PI * K(b(j)) * (dMw(b(j)) * Rho(b(j),1) * x2(b(j),1) ...
                                 + dMnw(b(j)) * Rho(b(j),2) * x2(b(j),2)) * (Prod(i).p - p(b(j)));
            end
        case('RateConstrained')
            for j=1:length(b)
                J1S(b(j),b(j)) = J1S(b(j),b(j));
                J2S(b(j),b(j)) = J2S(b(j),b(j));
            end
    end
end

end