%% Homework 4

R_saturn  = 60268;
Gm_saturn = 3.794*10^7;

r1 = [-720000; 670000; 310000];
v1 = [2.160; -3.360; 0.620];

% norm of the vectors
r1_norm = norm(r1);
v1_norm = norm(v1);

% angular momentum
h1      = cross(r1, v1);
h1_norm = norm(h1);

% eccentricity vector
e1      = cross(v1, h1) / Gm_saturn - r1 / r1_norm;
e1_norm = norm(e1)

% true anomoly st impact - radius of saturn plugged in!
trueAnomoly2 = acosd(h1_norm^2/Gm_saturn/(e1_norm*R_saturn) - 1/e1_norm )

if dot(r1,v1) > 0
    disp('true anomoly positive')
else
    disp('nrgative true anomoly')
    trueAnomoly2 = -trueAnomoly2
end


% find starting true anomaly - position of starting
trueAnomoly1 = acosd(h1_norm^2/Gm_saturn/(e1_norm*r1_norm) - 1/e1_norm )

if dot(r1,v1) > 0
    disp('true anomoly positive')
else
    disp('nrgative true anomoly')
    trueAnomoly1 = -trueAnomoly1
end

% All above work is from HW3

% obtain f and g funcitons

% change in true anomaly
deltaTA = trueAnomoly2 - trueAnomoly1;

p = h1_norm^2/Gm_saturn;

f = 1 - r1_norm/p * (1-cosd(deltaTA));

g = r1_norm*R_saturn/sqrt(Gm_saturn*p) * sind(deltaTA);

f_dot = sqrt(Gm_saturn/p) * tand(deltaTA/2) * ((1-cosd(deltaTA))/p - 1/R_saturn - 1/r1_norm);

g_dot = 1 - (r1_norm/p) * (1-cosd(deltaTA));

% now see the results of the f and g functions

r2 = f*r1 + g*v1

v2 = f_dot*r1 + g_dot*v1


% part b - want to find the time between each true anomlay

E1 = acosd((e1_norm + cosd(trueAnomoly1)) / (1 + e1_norm*cosd(trueAnomoly1)));
E1 = -E1;
% must be in radians!
E1 = deg2rad(E1)

E2 = acosd((e1_norm + cosd(trueAnomoly2)) / (1 + e1_norm*cosd(trueAnomoly2)));
E2 = -E2;
% must be in radians!
E2 = deg2rad(E2);

% find mean motion 
energy1 = v1_norm^2 / 2 - Gm_saturn / r1_norm;

a1 = -Gm_saturn / (2*energy1);

n = sqrt(Gm_saturn/(a1^3)); 

% Mean anomaly 

M1 = E1 - e1_norm*sin(E1);
M2 = E2 - e1_norm*sin(E2);

% time past periapsis
tp1 = M1 / n;
tp2 = M2 / n;

elaspedTime = tp2 - tp1