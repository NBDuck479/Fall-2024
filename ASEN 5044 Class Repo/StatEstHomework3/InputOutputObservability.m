function [x0] = InputOutputObservability(Udata, Ydata, F, G, H, M)

Udata = Udata'; 
Ydata = Ydata';

% observabilityMatrix = [H; H*F; H*F^2; H*F^3]

% Y{1} = Ydata(:,1)
% 
% Y{2} = Ydata(:,2) - H*G*Udata(:,1)
% 
% Y{3} = Ydata(:,3) - H*F*G*Udata(:,1) - H*G*Udata(:,2)
% 
% Y{4} = Ydata(:,4) - H*F^2*G*Udata(:,1) - H*F*G*Udata(:,2) - H*G*Udata(:,3)


% Loop to calculate Y{2} to Y{length(Ydata)}
for k = 2:length(Ydata)
    % Initialize the current value
    currentY = Ydata(:, k);
    
    % Subtract contributions from Udata
    for j = 1:k-1
        currentY = currentY - H * F^(k-1-j) * G * Udata(:, j);
    end
    
    % Also set Obsevability Matrix
    observabilityMatrix{k} = H*F^(k-1);
    
    % Store the result
    Y{k} = currentY;
end

% started with Y{2} in the algorithm so set Y{1} manually 
% Y{1} = Ydata(:,1);


% cell to mat
Y = cell2mat(Y');

observabilityMatrix = cell2mat(observabilityMatrix');


x0 = (observabilityMatrix'*observabilityMatrix)^(-1)*observabilityMatrix'*Y;
