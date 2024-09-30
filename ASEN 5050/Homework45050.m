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
trueAnomoly2 = acosd(h1_norm^2 / (e1_norm*R_saturn*Gm_saturn) - 1/e1_norm)

if dot(r1,v1) > 0
    disp('true anomoly positive')
else
    disp('nrgative true anomoly')
    trueAnomoly2 = -trueAnomoly2
end


% find starting true anomaly - position of starting
trueAnomoly1 = acosd(h1_norm^2/ (e1_norm*r1_norm*Gm_saturn) - 1/e1_norm)

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

f = 1 - (r1_norm/p) * (1-cosd(deltaTA));

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

%% part 2 

r1 = [5.352950*10^6; 7.053778*10^5; -4.059700*10^5]; 
v1 = [-4.164248; 1.963690; 3.191257*10^-1]; 

Gm_Jupiter = 1.268*10^8;
R_Jupiter  = 71492; 

% part a - true anomaly and eccen anomaly 

% norm of the vectors
r1_norm = norm(r1);
v1_norm = norm(v1);

% angular momentum
h1      = cross(r1, v1);
h1_norm = norm(h1);

% eccentricity vector
e1      = cross(v1, h1) / Gm_Jupiter - r1 / r1_norm;
e1_norm = norm(e1);

% true anomoly st impact - radius of saturn plugged in!
trueAnomoly1 = acosd(h1_norm^2 / (e1_norm*r1_norm*Gm_Jupiter) - 1/e1_norm);

if dot(r1,v1) > 0
    disp('true anomoly positive')
else
    disp('nrgative true anomoly')
    trueAnomoly1 = -trueAnomoly1
end

% eccentric anomaly
E1 = acosd((e1_norm + cosd(trueAnomoly1)) / (1 + e1_norm*cosd(trueAnomoly1)));
E1 = -E1;
% must be in radians!
E1 = deg2rad(E1)



% part b

% find mean motion 
energy1 = v1_norm^2 / 2 - Gm_Jupiter / r1_norm;

a1 = -Gm_Jupiter / (2*energy1);

meanMotion = sqrt(Gm_Jupiter/(a1^3)); 

% line of nodes
LON = cross([0,0,1], h1); 

% fid the arg of periapsis 

w1 = acosd(dot(LON,e1) / (norm(LON)*e1_norm))

if dot(e1,[0,0,1]) > 0
    disp('positive w')
else
    disp('negative w')
    w1 = -w1;
end 

% trueAnomaly
trueAnomaly2 = 180 - w1

% eccentric anomaly
E2 = acosd((e1_norm + cosd(trueAnomaly2)) / (1 + e1_norm*cosd(trueAnomaly2)));

% must be in radians!
E2 = deg2rad(E2)


% part c

% Mean anomaly 

M1 = E1 - e1_norm*sin(E1);
M2 = E2 - e1_norm*sin(E2);

% time past periapsis
tp1 = M1 / meanMotion;
tp2 = M2 / meanMotion;

elaspedTime = tp2 - tp1;

% need to find r2 norm for F and G functions
r2_norm = h1_norm^2/Gm_Jupiter / (1+e1_norm*cosd(trueAnomaly2));

% part d 

deltaTA = trueAnomaly2 - trueAnomoly1;

p = h1_norm^2/Gm_Jupiter;

f = 1 - r1_norm/p * (1-cosd(deltaTA));

g = r1_norm*r2_norm/sqrt(Gm_Jupiter*p) * sind(deltaTA);

f_dot = sqrt(Gm_Jupiter/p) * tand(deltaTA/2) * ((1-cosd(deltaTA))/p - 1/r2_norm - 1/r1_norm);

g_dot = 1 - (r1_norm/p) * (1-cosd(deltaTA));

r2 = f*r1 + g*v1

v2 = f_dot*r1 + g_dot*v1



% part e 

SMA = a1;
ECCEN = e1_norm;
MU = Gm_Jupiter;
tp = 1000; 
E = KeplersEQsolver(tp, SMA, ECCEN, MU)