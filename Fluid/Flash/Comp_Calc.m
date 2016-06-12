function [x,y] = Comp_Calc(P,T,b,Tb,z,Tol)

%Takes pressure as Pa
%Takes temperature as K
%b in K. this is the slope of the straight line connecting the critical point and the atmospheric boiling point on a log vapor pressure vs liT plot 
%Takes Tb as K
%Z is unitless

%This solves for the mole fractions in each phase for a two phase multi
%component system. This uses the Standing (1979) correlation and is really
%only valid for pressure less than 1000 psia and temperatures less than 200
%deg F

%Barnaby Fryer

%% Inputs that will later be part of the function

fvguess = .5;           %Initial guess for vapor fraction
tol = Tol;              %Defines the tolerance
alpha = .1;             %Correction factor to help convergence (higher convergences faster but has more chance to fail)
hurry = 100;
breakers = 0;

%% Basic Calculations
T = T*(9/5);            %From K to R
b = b*(9/5);            %From K to R
Tb = Tb*(9/5);          %From K to R
P = P*0.000145037738;   %From Pa to psi

%% Finding K

F = b .* ((1./Tb) - (1/T));                     %Finds F factor as per Standing 1979
a = 1.2 + 4.5 * 10^-4 * P + 15 * 10^-8 * P^2;   %a coefficient
c = 0.89 - 1.7 * 10^-4 * P - 3.5 * 10^-8 * P^2; %c coefficient
k = (10.^(a + c * F))/P;                        %K as per Standing 1979

%% Checking if in two phase region

BubCheck = z .* k;                  %Checks if below bubble point
BubCheck = sum(BubCheck);

if BubCheck < 1
    %display('Below bubble point. Only liquid present. Need higher temperature of lower pressure.')
    x = z;
    y = zeros(1,length(z));
    y(1,2) = 1;
    breakers = 1;
end

DewCheck = z ./ k;                  %Checks if above dew point
DewCheck = sum(DewCheck);

if DewCheck < 1
    %display('Above dew point. Only vapor present. Need lower temperature of higher pressure.')
    y = z;
    x = zeros(1,length(z));
    x(1,1) = 1;
    breakers = 1;
end

%% Solving for fv

fv = fvguess;                           %Gives the initial guess to the function
error = 1;
counter = 0;

while error > tol
    hi = (z .* (k)) ./ (fv .* (k - 1) + 1);             %Finds hi for each component
    dhi = (z .* (k - 1).^2) ./ ((fv .* (k - 1) + 1).^2);%Finds the derivative of hi for each component
    h = sum(hi) - 1;                                    %Sums hi (should be zero)
    dh = -sum(dhi);                                     %Sums derivative (used to correct fv)
    
    fvnew = alpha * (-h / dh) + fv;                     %Updates fv for error calculation
    error = abs((fvnew - fv) / fv);                     %Finds error associated with this result
    fv = fvnew;                                         %Updates fv for next iteration
    if fv < 0                                           %Prevents function from crashing or yielding unrealistic values (fv must be between 0 and 1)
        fv = rand();
    end
    if fv > 1                                           %Prevents function from crashing or yielding unrealistic values (fv must be between 0 and 1)
        fv = rand();
    end
    counter = counter + 1;
    if counter > hurry
        alpha = alpha / 2;
        counter = 5;
    end
end

%% Solve for x's and y's

if breakers ~= 1
    x = (z) ./ (fv .* (k - 1) + 1);         %Solves for mole fractions in liquid phase
    x(1,1) = 1 - x(1,2);
    y = k .* x;                             %Solves for mole fractions in gas phase
    y(1,2) = 1 - y(1,1);
end

end

