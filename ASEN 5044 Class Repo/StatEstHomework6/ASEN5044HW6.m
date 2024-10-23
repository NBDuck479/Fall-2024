%% Stat Est HW6

% part a in writing


%% part b 
t = 0; 
deltaT = 0.2;
%PSD
W = 10; 

% Matrices
A = [0 1; -100 -10];
Gamma = [0; 1]; 
% maybe this istn actually B, but Gamma?

% F matrix is same as earlier in semester
F = expm(A*deltaT)

% find G 
A_hat = [A, Gamma; 0 0 0];

A_aug = expm(A_hat*deltaT); 

G = A_aug(1:2, 3);

% Van Loan method for solving for Q
Z_top = [-A Gamma*W*Gamma'];
Z_bot = [zeros(2), A'];

Z = [Z_top; Z_bot];

% matrix exp of Z
z_expm = expm(Z*deltaT);

% grab the top right of matrix
FinvQ = z_expm(1:2, 3:4); 

Ftrans = z_expm(3:4, 3:4); 

% solve for Q
Q = Ftrans' * FinvQ;

%% Part d 