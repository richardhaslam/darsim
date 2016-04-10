%Build FIM Jacobian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2016
%Last modified: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Jop, Jwp, Jos, Jws] = AddWellsToJacobian(Jop, Jwp, Jos, Jws, Inj, Prod, K, p, Mw, Mo, dMw, dMo)
%Injectors
for i=1:length(Inj)
    a = Inj(i).cells;
    for j=1:length(a)
        Jop(a(j),a(j)) = Jop(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mo;
        Jwp(a(j),a(j)) = Jwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mw;
    end
end
%Producers
for i=1:length(Prod)
    b = Prod(i).cells;
    for j=1:length(b)
        Jop(b(j),b(j)) = Jop(b(j),b(j)) + Prod(i).PI*K(b(j)).*Mo(b(j));
        Jwp(b(j),b(j)) = Jwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*Mw(b(j));
        Jos(b(j),b(j)) = Jos(b(j),b(j)) - Prod(i).PI*K(b(j)).*(Prod(i).p-p(b(j))).*dMo(b(j));
        Jws(b(j),b(j)) = Jws(b(j),b(j)) - Prod(i).PI*K(b(j)).*(Prod(i).p-p(b(j))).*dMw(b(j));
    end
end

end