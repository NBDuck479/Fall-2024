%% HW2

% Problem 1
Gm_mars = 4.305*10^4; 

r = [3.62067*10^3; -3.19925*10^2; -4.20645*10^2]; 
v = [-4.28843*10^-1; -3.00176*10^-2; -3.39801]; 

r_norm = norm(r);
v_norm = norm(v);

h      = cross(r,v);
h_norm = norm(h);

% find a,e,i,RAAN,w,True Anomoly

% calulate the specific energy
energy = v_norm^2 / 2 - Gm_mars / r_norm

a = -Gm_mars / (2*energy)

e_norm = sqrt(1 + (2*h_norm^2*energy)/(Gm_mars^2))

% inclination
i = acosd(h(3)/h_norm)

% line of nodes
n = cross([0,0,1], h)
n_norm = norm(n);

% RAAN
RAAN = acosd(dot(n,[1,0,0]) / n_norm)

if dot(n, [0,1,0]) > 0
    disp('positive RAAN')
else
    disp('negavtie RAAN')
end

% eccentricity vector
e = cross(v, h) / Gm_mars - r / r_norm;

% true anomoly
trueAnomoly = acosd(dot(r,e) / (r_norm * e_norm))

if dot(r,v) > 0 
    disp('true anomoly positive')
else
    disp('nrgative true anomoly')
    trueAnomoly = -trueAnomoly;
end

% arg of periapsis 
w = acosd(dot(n,e) / (n_norm*e_norm))

if dot(e,[0,0,1]) > 0
    disp('positive w')
else
    disp('negative w')
    w = -w;
end

theta = trueAnomoly + w;


C1 = [cosd(RAAN)*cosd(theta)-sind(RAAN)*cosd(i)*sind(theta) -cosd(RAAN)*sind(theta)-sind(RAAN)*cosd(i)*cosd(theta) sind(RAAN)*sind(i)];
C2 = [sind(RAAN)*cosd(theta)+cosd(RAAN)*cosd(i)*sind(theta) -sind(RAAN)*sind(theta)+cosd(RAAN)*cosd(i)*cosd(theta) -cosd(RAAN)*sind(i)];
C3 = [sind(i)*sind(theta) sind(i)*cosd(theta) cosd(i)];

C = [C1; C2; C3]


r_rth = C' * r

v_rth = C' * v


theta = -65.16;
RAAN = 175.0805;
i = 91.1242;

C1 = [cosd(RAAN)*cosd(theta)-sind(RAAN)*cosd(i)*sind(theta) -cosd(RAAN)*sind(theta)-sind(RAAN)*cosd(i)*cosd(theta) sind(RAAN)*sind(i)];
C2 = [sind(RAAN)*cosd(theta)+cosd(RAAN)*cosd(i)*sind(theta) -sind(RAAN)*sind(theta)+cosd(RAAN)*cosd(i)*cosd(theta) -cosd(RAAN)*sind(i)];
C3 = [sind(i)*sind(theta) sind(i)*cosd(theta) cosd(i)];



 rp = h_norm^2/Gm_mars / (1+e_norm*cosd(0))

 
 vr = Gm_mars/h_norm*e_norm*sind(0)

 vt = Gm_mars/h_norm*(1+e_norm*cosd(0))
 
 r_XYZ = C * [3613.69; 0; 0]
v_XYZ = C * [0; 3.4679; 0]
 
 %% Problem 2
clear all; clc
   
% a) 

n_hat = [0.6428; -0.7660; 0];
h_hat = [-0.3237; -0.2717; 0.9063];
e     = [0.0475; 0.3755; 0.1295];
e_hat = e / norm(e);


acosd(dot(n_hat,e) / (norm(n_hat)*norm(e)))




n_norm = norm(n_hat)
e_norm = norm(e)

i = acosd(h_hat(3)/norm(h_hat))


RAAN = acosd(dot(n_hat,[1,0,0]) / norm(n_hat))

if dot(n_hat, [0,1,0]) > 0
    disp('positive RAAN')
else
    disp('negavtie RAAN')
    RAAN = -RAAN; 
end

% arg of periapsis 
w = acosd(dot([0.6428; -0.7660; 0], [0.0475; 0.3755; 0.1295]) / (norm([0.6428; -0.7660; 0]) * norm([0.0475; 0.3755; 0.1295])))

if dot(e,[0,0,1]) > 0
    disp('positive w')
else
    disp('negative w')
end


% part b) 
r = 4070.6;

r_n = r * n_hat
r_n_norm = norm(r_n);
r_n_hat = r_n / r_n_norm;

% true anomoly
trueAnomoly = asind(dot(r_n_hat,[0,0,1]) / sind(i)) - w


SMA = r_n_norm * (1 + e_norm*cosd(trueAnomoly)) / (1-e_norm^2)
