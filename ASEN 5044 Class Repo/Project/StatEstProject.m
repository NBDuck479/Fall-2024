%% 5044 Project Stat OD start 
% I'll try to use camel case for variables 
% Use spacing and align things for ease of reading 


%% Given constants and parameters: 
% Earth Gravitational Parameter
muEarth = 398600; 

% Radius of Earth
earthRadius = 6378; % km

% Rotation of Earth 
omegaEarth = 2*pi / 86400; % rad/sec



%% *Part I* 

% initial states
xInitial = 6678; 
yInitial = 0;

% intial position
r0 = sqrt(xInitial^2 + yInitial^2);

% initial velocity
xDotInitial = 0; 
yDotInitial = r0*sqrt(muEarth/(r0^3)); 


% Nominal Jacobian - I'm evaluating it at the given initial conditions (I guess that's what I do?)
nomJacobian = [0, 1, 0, 0; ...
 -muEarth/(r0^3), 0, 0, 0; ...
               0, 0, 0, 1]


