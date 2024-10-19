%% *Homework 5 script*

mu_mars = 4.305*10^4; 
mu_moon = 4902.799; 
mu_sun  = 1.32712428*10^11; 
mu_sat  = 3.794*10^7; 
R_mars  = 3397.2; 
R_moon  = 1738; 


%% problem 1

% part - a 

SMA1 = 7045; 
e1   = 0.23; 
TA   = -142; 

% energy of orbit
energy = -mu_moon  / (2*SMA1);

% specific angular momentum 
h = sqrt(mu_moon*SMA1*(1-e1^2));

% radial velcoty
Vr = mu_moon/h*(e1*sind(TA));

%theta velo
Vt = mu_moon/h*(1+e1*cosd(TA));

% velcoity
v1_rth = [Vr; Vt; 0];

% part - c
 deltaV_rth = [0.3; -.1; 0];
 v2 = v1_rth + deltaV_rth;
 
 
 % part - d
 % get a2 e2 TA2 
 
 % conic for position
 r1 = h^2/mu_moon / (1+e1*cosd(TA))
 
 % get energy - position stays the same at manuever 
 energy2 = norm(v2)^2/2 - mu_moon/r1;
 
 % SMA2 
 SMA2 = -mu_moon / (2*energy2);
 
 % angular momentun 
 h2 = cross([r1,0,0], v2)
 
 % ecce
 e2 = sqrt(1 + (2*norm(h2)^2*energy2) / (mu_moon^2))
 
 % TA2
 TA2 = acosd(norm(h2)^2/(e2*r1*mu_moon) - 1/e2)
 
 if dot([r1;0;0], v2 ) > 0 
    disp('true anomoly positive')
else
    disp('nrgative true anomoly')
    TA2 = -TA2;
 end

 
 % part e - change in argument of periapsis!
 
deltaW = TA - TA2;
 % part f - mass 

 Mass      = 1224 % kg
 deltaV_ms = norm(deltaV_rth) * 1000 % m/s
 g         = 9.81; %m/s^2
 Isp       = 212;
 
 mp = Mass*(1 - exp(-deltaV_ms/(Isp*g)))
 
 %% problem 2
 
 r1  = 6500; 
 E1  = pi/2; 
 rp1 = 5915;
 
 rp2 = 5712; 
 ra2 = 7888; 
 
 % TA 
 LHS = (6500/5915)^2;
 
 syms ecc 
 eqn = (LHS+1)*ecc^2 - 2*LHS*ecc + LHS-1 ==0; 
 
 e1 = solve(eqn, ecc);
 e1 = 1719/18281;
 
 % spec angular momentum 
 h1 = sqrt(mu_mars*rp1*(1+e1*cosd(0))); 
 
 % True Anomaly 1
 TA1 = acosd(h1^2/(mu_mars*r1*e1) - 1/e1)
 
 % Vr
  V1r = (mu_mars/h1)*e1*sind(TA1);
  V1t = (mu_mars/h1)*(1+e1*cosd(TA1)); 
  
 % Stacked Vrth vector
 V1rth = [V1r; V1t; 0]
  
 
 
 % part b
 
 SMA2 = .5*(ra2 + rp2);
 e2   = (ra2 - rp2) / (ra2 + rp2);
 
 % specific angular momentum
 h2 = sqrt(mu_mars*SMA2*(1-e2^2));
 
 % TA2
 TA2 = acosd(h2^2/(mu_mars*r1*e2) - 1/e2);
 
 % velocity
 V2r = mu_mars/h2*e2*sind(TA2);
 V2t = mu_mars/h2*(1+e2*cosd(TA2));
 
 V2rth = [V2r; V2t; 0]
 
 % change in velocity answer
 deltaV = V2rth - V1rth
 
 norm(deltaV)
 
 
 %% problem 3
 
 AU         = 149597870.7; % km
 SMA_earth  = AU * 1.0000010178;
 SMA_saturn = AU * 9.554909595;
 
 
% velocity of s/c in earth orbit
V1 = sqrt(mu_sun/SMA_earth)

% transfer arc
SMA_trans = .5*(SMA_earth + SMA_saturn);

% energy of transfer
energy_trans = -mu_sun/(2*SMA_trans); 

% velocity at start of transfer arc
V1_trans = sqrt(2*(energy_trans + mu_sun/SMA_earth));

% velocity at end of transfer arc 
V2_trans = sqrt(2*(energy_trans + mu_sun/SMA_saturn));

% Velocity needed for saturn orbit
V2 = sqrt(mu_sun / SMA_saturn);

% get deltaV for each manuever
deltaV1_trans = V1_trans - V1;
deltaV2_trans = V2 - V2_trans;

total_deltaV = deltaV1_trans + deltaV2_trans

% time of flight
TOF = pi*sqrt(SMA_trans^3 / mu_sun)



% part - b
n2 = sqrt(mu_sun / (SMA_saturn^3)); 

alpha = n2 * TOF;

psi = pi - alpha;


% part - c
RB = AU * 11; %km

% init veloc
Vi1 = sqrt(mu_sun/SMA_earth);

% first transfer velo needed
SMA_trans1 = .5*(SMA_earth + RB); 

energy_trans1 = -mu_sun/(2*SMA_trans1); 

V1f = sqrt(2*(energy_trans1 + mu_sun/SMA_earth));

% calculate the first delta V
deltaV1 = V1f - Vi1;

% Velocity at RB of transfer orbit V2i
V2i = sqrt(2*(energy_trans1 + mu_sun/RB))

% velocity needed for second transfer 
SMA_trans2 = 0.5*(RB + SMA_saturn);

% energy of second transfer orbit
energy_trans2 = -mu_sun/(2*SMA_trans2); 

% velocity for second transfer 
V2f = sqrt(2*(energy_trans2 + mu_sun/RB));

% deltaV needed for the second transfer
deltaV2 = V2f - V2i;


% finally, the last impulse needed to complete bi-ellitpc transfer
V3i = sqrt(2*(energy_trans2 + mu_sun/SMA_saturn)); 

% energy of saturns orbit 
energy_saturn = -mu_sun/(2*SMA_saturn);

% velocity for saturns orbit
V3f = sqrt(mu_sun/SMA_saturn);

% deltaV needed for last impulse
deltaV3 = V3f - V3i;


% total deltaV
deltaVtot = abs(deltaV1) + abs(deltaV2) + abs(deltaV3)

% time of flight for each transfer to get total
TOF1 = pi*sqrt(SMA_trans1^3 / mu_sun);

TOF2 = pi*sqrt(SMA_trans2^3 / mu_sun);

totalTOF = TOF1+TOF2

% part d 

