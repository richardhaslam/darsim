function [Status, SinglePhase] = SimpleFlashCalculator(Status)
% Define SinglePhase objects
SinglePhase = zeros(length(Status.p), 1);
z = Status.z;
x = Status.x;
ni = Status.ni;
N = length(ni);
nc = 3;

k(:, 1) = ones(N, 1) * 1.5;
k(:, 2) = ones(N, 1) * 0.5;
k(:, 3) = ones(N, 1) * 0.1;

%% 2 Chek if we are in 2 phase region
% 2.a: checking if it 's all liquid: checks if mix is below bubble
BubCheck = z .* k;
BubCheck = sum(BubCheck, 2);
SinglePhase(BubCheck - 1 <= 1e-6) = 2; % It s all liquid

% 2.b: checking if it 's all vapor: checks if mix is above dew
% point
DewCheck = z ./ k;
DewCheck = sum(DewCheck, 2);
SinglePhase(DewCheck - 1 < 1e-6) = 1;  % It s all vapour

% Solve with Newton
converged = 0;
for i=1:N
    switch SinglePhase(i)
        case(0)
            while ~converged
                Residual = zeros(2*nc+1,1);
                J = zeros(2*nc+1);
                for c=1:nc
                    Residual(c) = x(i, (c-1)*2+1) - k(i, c) * x(i, (c-1)*2+2); 
                    Residual(nc+c) = z(i, c) - ni(i)*x(i, (c-1)*2+1) -(1-ni(i))*x(i, (c-1)*2+2);                 
                    Residual(2*nc+1) = Residual(2*nc+1) + x(i, (c-1)*2+1) - x(i, (c-1)*2+2);
                    % x_cv - k_c x_cl = 0
                    J (c, (c-1)*2+1:(c-1)*2+2) = [1, -k(i,c)];
                    % z_c - ni*x_cv - (1-ni)*x_cl = 0
                    J (nc + c, (c-1)*2+1:(c-1)*2+2) = [-ni(i), ni(i)-1];
                    J (nc + c, end) = x(i, (c-1)*2+2) - x(i, (c-1)*2+1); 
                    % Sum x_cv-x_cl = 0
                    J (2*nc + 1, (c-1)*2+1:(c-1)*2+2) = [1, -1];
                end
                dx = -J\Residual;
                x(i, :) = x(i, :) + dx(1:end-1)';
                ni(i) = ni(i) + dx(2*nc+1);
                if norm(Residual, inf) < 1e-6
                    converged = 1;
                end
            end
        case(1)
            % It s all vapour
            for c=1:nc-1
                x(i,(c-1)*2+1) = z(i, c);
                x(i,(c-1)*2+2) = 1; % This is for Stability
            end
            ni(i) = 1;
            converged = 1;
        case(2)
            % It s all liquid
            for c = 1:nc
                x(i, (c-1)*2+2) = z(i,c);
                x(i, (c-1)*2+1) = 1; % This is for Stability
            end
            ni(i) = 0;
            converged =1;
    end
    if ~converged
        disp('Warning: Flash did not fully converge!')
    end
end
% Copy it into Status object
Status.x = x;
Status.ni = ni;
end