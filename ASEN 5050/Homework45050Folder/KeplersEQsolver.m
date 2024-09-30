function [E] = KeplersEQsolver(tp, SMA, ECCEN, MU)
% This functio will iteratively solve keplers equation to find the
% eccentric anomaly 

% Define initial guess for solver 
% Common Approach of guessing mean anomaly is eccentric anomaly 

meanMotion = sqrt(MU / (SMA^3)); 

meanAnomaly = meanMotion*tp; 

% set initial condition 
E0 = deg2rad(meanAnomaly); 

% Keplers Eq'n set to zero

g = @(E) E - ECCEN * sin(E) - meanMotion; 

g_prime = @(E) 1 - ECCEN * cos(E); 

% Set max iteration limit
count = 1;

while g(E0) > 10^(-10) && count < 15

    E_future = E0 - g(E0)/g_prime(E0);
    
    E0 = E_future;
    
    count = count + 1;
    
end

% give the answer

E = E0; 