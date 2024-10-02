function [R_XYZ, V_XYZ] = Orbital2Cartesian(SMA, ecc, i, RAAN, w, TA, mu)
% MOST IMPORTANT THAT THE ORBITAL ELEMENTS SHOULD BE -180 to 180 DEGREES!!

% this function computes orbital elements to cartesian for specified mu. 
% all the units you'd expect 

% specific energy
energy = -mu / (2*SMA); 

% specific angular momentum
h = sqrt((ecc^2-1)*mu^2 / (2*energy)); 

% velocity in rth frame: 

Vr = mu/h * ecc * sind(TA); 

Vt = mu/h * (1 + ecc*cosd(TA));

% velocity in rth frame 
V_rth = [Vr; Vt; 0]; 

% position in the rth frame 
r = h^2/mu / (1+ecc*cosd(TA));

% position in the rth frame
R_rth = [r; 0; 0]; 

% Now transform the frame to the primary body inertial frame 

theta = TA + w; 

C1 = [cosd(RAAN)*cosd(theta)-sind(RAAN)*cosd(i)*sind(theta) -cosd(RAAN)*sind(theta)-sind(RAAN)*cosd(i)*cosd(theta) sind(RAAN)*sind(i)];
C2 = [sind(RAAN)*cosd(theta)+cosd(RAAN)*cosd(i)*sind(theta) -sind(RAAN)*sind(theta)+cosd(RAAN)*cosd(i)*cosd(theta) -cosd(RAAN)*sind(i)];
C3 = [sind(i)*sind(theta) sind(i)*cosd(theta) cosd(i)];

C = [C1; C2; C3]; 

% Now convert the rth to the inertal frame

R_XYZ = C * R_rth; 
V_XYZ = C * V_rth;