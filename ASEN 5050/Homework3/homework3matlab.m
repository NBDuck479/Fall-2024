%% Homework 3 

r1 = [-720000; 670000; 310000]; 
v1 = [2.160; -3.360; 0.620]; 

Gm_saturn = 3.794*10^7; 
R_saturn  = 60268; 

% norm of the vectors
r1_norm = norm(r1);
v1_norm = norm(v1);

% get position and velcity at impact 

% angular momentum 
h1      = cross(r1, v1); 
h1_norm = norm(h1);

% inclination
i = acosd(h1(3)/h1_norm)

% line of nodes
n1 = cross([0,0,1], h1);
n1_norm = norm(n1);

% RAAN
RAAN1 = acosd(dot(n1,[1,0,0]) / n1_norm)

if dot(n1, [0,1,0]) > 0
    disp('positive RAAN')
else
    disp('negavtie RAAN')
    RAAN1 = -RAAN1;
end

% eccentricity vector
e1      = cross(v1, h1) / Gm_saturn - r1 / r1_norm;
e1_norm = norm(e1)

% true anomoly
trueAnomoly1 = acosd(h1_norm^2/Gm_saturn/(e1_norm*R_saturn) - 1/e1_norm )

if dot(r1,v1) > 0 
    disp('true anomoly positive')
else
    disp('nrgative true anomoly')
    trueAnomoly1 = -trueAnomoly1;
end

% arg of periapsis 
w1 = acosd(dot(n1,e1) / (n1_norm*e1_norm))

if dot(e1,[0,0,1]) > 0
    disp('positive w')
else
    disp('negative w')
    w1 = -w1;
end

% calulate the specific energy
energy1 = v1_norm^2 / 2 - Gm_saturn / r1_norm

a1 = -Gm_saturn / (2*energy1)

% now that all the orbital elements are found, we can say that impact will
% occur at rp and calculate everything else! 

% We can get velocity in the r theta h frame! 

theta_star = trueAnomoly1;

Vr = Gm_saturn/h1_norm * e1_norm * sind(theta_star); 

Vtheta = Gm_saturn/h1_norm * (1 + e1_norm*cosd(theta_star)); 

V_rth = [Vr; Vtheta; 0];

% position in rth is rp all in r direction 
% HAVING rp BEING CALCULATED BUT COULD PROB ALSO HAVE IT SATURN RADIUS 
rp = R_saturn;

R_rth = [rp; 0; 0];

% Now get the position and velocity from rth to XYZ coordinates 
theta = trueAnomoly1 + w1; 

C1 = [cosd(RAAN1)*cosd(theta)-sind(RAAN1)*cosd(i)*sind(theta) -cosd(RAAN1)*sind(theta)-sind(RAAN1)*cosd(i)*cosd(theta) sind(RAAN1)*sind(i)];
C2 = [sind(RAAN1)*cosd(theta)+cosd(RAAN1)*cosd(i)*sind(theta) -sind(RAAN1)*sind(theta)+cosd(RAAN1)*cosd(i)*cosd(theta) -cosd(RAAN1)*sind(i)];
C3 = [sind(i)*sind(theta) sind(i)*cosd(theta) cosd(i)];

C = [C1; C2; C3]; 

% C should be able to spit out XYZ coordinates 

% Dont forget to discuss the results! 




%% Problem 2: 
% 
mu_mars = 42828.314258067; 


SMA  = 6463.8; 
ecc  = 0.45454; 
i    = 74.924; 
RAAN = 1.2410; 
w    = 353.31 - 360; 
TA   = 199.38 - 360; 

% orbital Period 

period = 2*pi*sqrt(SMA^3/mu_mars); % seconds
period = period/60/60; % Period is in hours! 

rp = SMA*(1 - ecc); 
% SHOULD I SUBTRACT RADIUS OF MARS? FOR ALTITUDE PART? 


% Part c
[R_XYZ, V_XYZ] = Orbital2Cartesian(SMA, ecc, i, RAAN, w, TA, mu_mars)

% part G 

ECCMars = dlmread('ECCmars');

for i = 1:length(ECCMars)
    ECCmars_norm(i) = norm(ECCMars(:,2));
end

fig = 1;
figure(fig)
plot(0.00289403075109033:0.00289403075109033:1.826133403938, ECCmars_norm)
grid on
ylabel('Eccentricity')
xlabel('Time [Days]')
title('Mars Point Mass - 10 orbits')
fig = fig + 1; 

HMAGmars = dlmread('HMAGmars');

for i = 1:length(HMAGmars)
    HMAGmars_norm(i) = norm(HMAGmars(:,2));
end

figure(fig)
plot(0.00289403075109033:0.00289403075109033:1.826133403938, HMAGmars_norm)
grid on 
ylabel('Specific Angular Momentum')
xlabel('Time [Days]')
title('Mars Point Mass - 10 orbits')
fig = fig+1; 

% part j 
ECCmarsHF = dlmread('ECCmarsHF');


for i = 1:length(ECCmarsHF)
    ECCmarsHF_norm(i) = norm(ECCmarsHF(i,2));
end

figure(fig)
plot(linspace(0.000694444443070097,18.2613340393473, length(ECCmarsHF_norm)), ECCmarsHF_norm)
grid on
ylabel('Eccentricity')
xlabel('Time [Days]')
title('Mars High Fidelity - 100 orbits')
fig = fig + 1; 

HMAGmarsHF = dlmread('HMAGmarsHF');

for i = 1:length(HMAGmarsHF)
    HMAGmarsHF_norm(i) = norm(HMAGmarsHF(i,2));
end

figure(fig)
plot(linspace(0.000694444443070097,18.2613340393473, length(HMAGmarsHF_norm)), HMAGmarsHF_norm)
grid on
ylabel('Specific Angular Momentum')
xlabel('Time [Days]')
title('Mars High Fidelity - 100 orbits')