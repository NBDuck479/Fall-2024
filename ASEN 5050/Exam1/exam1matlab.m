%% exam 1

mu_mars = 4.305*10^4; 
R_mars  = 3397.2; 

%% Probblem 1
 h_norm  = 2.4799*10^4; 
 r1_norm = 15000; 
 v1_norm = 1.8111; 
 r2_norm = 19000; 
 
 % part a - SMA 
 energy = v1_norm^2/2 - mu_mars/r1_norm;
 
 % use energy to solve for SMA
 
  SMA = -mu_mars / (2*energy);
  
  
  % part b - Eccentricity 
  
  ECCEN = sqrt(1 + (2*h_norm^2*energy)/(mu_mars^2));
  
  
  % part c - eccentric and true anomaly at t1
  
  TA1 = acosd(h_norm^2/(r1_norm*mu_mars*ECCEN) - 1/ECCEN);
  % Sign check on TA1 - negative fpa so negative TA1!
  TA1 = -TA1;
  
  % Eccentric anoamly - IN RADIANS!!!!
  E1 = acosd((ECCEN+cosd(TA1)) / (1+ECCEN*cosd(TA1)));
  E1 = -deg2rad(E1);
  
  
  
  % part d - TA2 and E2
  TA2 = acosd(h_norm^2/(r2_norm*mu_mars*ECCEN) - 1/ECCEN);
  % negative TA2 
  TA2 = -TA2;
  
  % E2
  E2 = acosd((ECCEN+cosd(TA2)) / (1+ECCEN*cosd(TA2)));
  E2 = -deg2rad(E2);
  
  
  
  % part e - Time differnece between both points 
  
  meanMotion    = sqrt(mu_mars/(SMA^3));
  meanAnonmaly1 = E1 - ECCEN*sin(E1);
  meanAnonmaly2 = E2 - ECCEN*sin(E2);
  
  % time until periapsis 
  tp1 = meanAnonmaly1 / meanMotion;
  tp2 = meanAnonmaly2 / meanMotion;
  
  % Calculate Period
  Period = 2*pi * sqrt(SMA^3 / mu_mars); 
  
  % interval diff
  intervalDiff = tp1 - tp2;
  
  % Now subtract from orbital period
  timeDiff = Period - intervalDiff;
  
  
  
  % part f - sketch variables 
  rp = SMA*(1-ECCEN);
  ra = SMA*(1+ECCEN);
  
  b = SMA*sqrt(1-ECCEN^2);
  
  semiLatus = h_norm^2/mu_mars;
  
  center = SMA*ECCEN;
  
  
  %% problem 2
  r3 = [-7.6650*10^3; 6.5468*10^3; -4.5740*10^2];
  v3 = [1.6334; 0.1226; -1.9455];
  
  r3_norm = norm(r3);
  v3_norm = norm(v3);
  
  
  % part a orbital elements 
  h      = cross(r3, v3);
  h_norm = norm(h);
  
  energy = v3_norm^2/2 - mu_mars/r3_norm;
  
  SMA = -mu_mars / (2*energy);
  
  ECCEN = sqrt(1 + (2*h_norm^2*energy)/(mu_mars^2));
  
  % eccentricity vector
    ECCEN_vec      = cross(v3, h) / mu_mars - r3 / r3_norm;
  
  TA3 = acosd(h_norm^2/(r3_norm*mu_mars*ECCEN) - 1/ECCEN);
  
  if dot(r3,v3) > 0 
    disp('true anomoly positive')
else
    disp('nrgative true anomoly')
    TA3 = -TA3;
  end
  
% inclination
i = acosd(h(3)/h_norm);

% line of nodes
LON = cross([0,0,1], h);
LON_norm = norm(LON);

% RAAN
RAAN = acosd(dot(LON,[1,0,0]) / LON_norm);

if dot(LON, [0,1,0]) > 0
    disp('positive RAAN')
else
    disp('negavtie RAAN')
    RAAN = -RAAN
end

% arg of periapsis 
w = acosd(dot(LON,ECCEN_vec) / (LON_norm*ECCEN));

if dot(ECCEN_vec,[0,0,1]) > 0
    disp('positive w')
else
    disp('negative w')
    w = -w
end



% part b  - 2 hours later 

meanMotion = sqrt(mu_mars/(SMA^3));

  E3 = acosd((ECCEN+cosd(TA3)) / (1+ECCEN*cosd(TA3)));
  E3 = -deg2rad(E3);
  
  meanAnonmaly3 = E3 - ECCEN*sin(E3);
  
  % time until periapsis 
  tp3 = meanAnonmaly3 / meanMotion;
  
  % time 4 past periapsis
  tp4 = tp3 + 2*60*60;
  
  % time 4 calcualtions 
  meanAnomaly4 = tp4*meanMotion;
  
  % use newtonraph to get E4 and then can get TA4 and keep going!
  
  [E4] = KeplersEQsolver(tp4, SMA, ECCEN, mu_mars);
  
  % true anomaly at 4
  TA4 = 2*atan(sqrt((1+ECCEN)/(1-ECCEN)) * tan(E4/2));
  TA4 = rad2deg(TA4)
  
  % now get position and velocity in the rth frame at t4! 
  Vr4 = mu_mars/h_norm * ECCEN * sind(TA4);
  Vt4 = mu_mars/h_norm * (1+ECCEN*cosd(TA4))
  
  % get position from conic equation: 
  r4 = h_norm^2/mu_mars / (1+ECCEN*cosd(TA4)); 
  
  V_rth = [Vr4; Vt4; 0]
  R_rth = [r4; 0; 0]
  
  % convert these to the inertial coordinate frame!
  
  theta = TA4 + w; 

C1 = [cosd(RAAN)*cosd(theta)-sind(RAAN)*cosd(i)*sind(theta) -cosd(RAAN)*sind(theta)-sind(RAAN)*cosd(i)*cosd(theta) sind(RAAN)*sind(i)];
C2 = [sind(RAAN)*cosd(theta)+cosd(RAAN)*cosd(i)*sind(theta) -sind(RAAN)*sind(theta)+cosd(RAAN)*cosd(i)*cosd(theta) -cosd(RAAN)*sind(i)];
C3 = [sind(i)*sind(theta) sind(i)*cosd(theta) cosd(i)];

C = [C1; C2; C3]; 

% now get the inertial coordinates! 

R_XYZ = C*R_rth
V_XYZ = C*V_rth

  

%% Problem 3
norm([-.64279, -0.76604, 0])

n = [0.64279, 0.76604, 0]

RAAN = acosd(dot(n,[1,0,0]) / (norm(n)))

if dot(n, [0,1,0]) > 0
    disp('positive RAAN')
else
    disp('negavtie RAAN')
    RAAN1 = -RAAN1;
end


w = acosd(dot(n,[0.02970, 0.97508, 0.21985]))

if dot([0.02970, 0.97508, 0.21985],[0,0,1]) > 0
    disp('positive w')
else
    disp('negative w')
    w = -w
end