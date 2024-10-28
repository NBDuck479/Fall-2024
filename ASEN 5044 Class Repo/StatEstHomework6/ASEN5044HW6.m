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

3/.2


%% Textbook problem 
R = [1 0; 0 4]; 
H = [1; 1]; 

inv(H'*inv(R)*H) * H' * inv(R)




%% Problem 3

% initial state
x0 = [1;1;1]; 

R = [8 5.15 6.5; 5.15 5 -4.07; 6.5 -4.07 50];

H = eye(3);

% resulting simulated noisy data 

for k = 1:100
    % Cholesky decomp for 'matrix square root'
    Sv = chol(R, 'lower');
    
    % random sample draw
    q = rand(3,1);
    
    y(1:3, k) = H*x0 + Sv*q;
    
end

fig = 1;
figure(fig)
plot(1:100, y(1, :)', 1:100, y(2, :)')
grid on 
axis equal
legend('yk(1)', 'yk(2)')
title('yk(1) and yk(2)')

fig = fig + 1; 

figure(fig)
plot(1:100, y(1, :)', 1:100, y(3, :)')
grid on 
axis equal
legend('yk(1)', 'yk(3)')
title('yk(1) and yk(3)')


fig = fig + 1; 

figure(fig)
plot(1:100, y(2, :)', 1:100, y(3, :)')
grid on
axis equal
legend('yk(2)', 'yk(3)')
title('yk(2) and yk(3)')


% part b - compute sample covariance matrix 
var(y(1, :))
var(y(2, :))
var(y(3, :))




% part - c
xhat = inv(H'*inv(R)*H) * H'*inv(R)