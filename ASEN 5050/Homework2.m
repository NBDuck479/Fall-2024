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
    RAAN = -RAAN; 
end

% eccentricity vector
e = cross(v, h) / Gm_mars - r / r_norm;

% true anomoly
trueAnomoly = acosd(dot(r,e) / (r_norm * e_norm))

if dot(r,v) > 0 
    disp('true anomoly positive')
else
    disp('nrgative true anomoly')
end

% arg of periapsis 
w = acosd(dot(n,e) / (n_norm*e_norm))

if dot(e,[0,0,1]) > 0
    disp('positive w')
else
    disp('negative w')
end