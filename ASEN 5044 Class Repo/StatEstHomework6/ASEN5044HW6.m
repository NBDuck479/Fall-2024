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


%% part b - compute sample covariance matrix 
var(y(1,:));
var(y(2,:));
var(y(3,:));
C = diag([var(y(1,:)), var(y(2,:)), var(y(3,:))]);

cov(y(1,:), y(2,:));
cov(y(2,:), y(3,:));
cov(y(1,:), y(3,:));


%% part - c

xhatLS = inv(H'*inv(R)*H) * H'*inv(R) * y(1:3,1)

% Do it with first 3 samples
H3 = [H; H; H];

R3 = blkdiag(R,R,R);

xhatLS3 = inv(H3'*inv(R3)*H3) * H3'*inv(R3) * [y(1:3,1); y(1:3,2); y(1:3,3)]

% Now with first 10 samples
H10 = [H; H; H; H; H; H; H; H; H; H];
R10 = blkdiag(R,R,R,R,R,R,R,R,R,R);

xhatLS10 = inv(H10'*inv(R10)*H10) * H10'*inv(R10) * [y(1:3,1); y(1:3,2); y(1:3,3); y(1:3,4); y(1:3,5); y(1:3,6); y(1:3,7); y(1:3,8); y(1:3,9); y(1:3, 10)]



%% Part - d

yK = importdata('hw6problem3data.csv'); 

Rd = blkdiag(R10, R10, R10);

Hd = [H10; H10; H10];

xhatLS = inv(Hd'*inv(Rd)*Hd) * Hd'*inv(Rd) * yK(:)

var(yK(1,:))
var(yK(2,:))
var(yK(3,:))

cov(yK(1,:), yK(2,:))
cov(yK(2,:), yK(3,:))
cov(yK(1,:), yK(3,:))

%% Part - e
xhatUWLS = inv(Hd'*Hd) * Hd' * yK(:)


%% part - f

% inital estimate
xhat_old = [0; 0; 0];

% intial error covariance
P_old = diag([100 100 100]);

% measurement something matrix
H = eye(3);

% recursive LLS
for i = 1:30
    
    Kk = P_old*H' * inv(H*P_old*H + R);
    
    xhat(:, i) = xhat_old + Kk*(yK(:,i) - H*xhat_old);
    
    P = (eye(3) - Kk*H) * P_old * (eye(3) - Kk*H)' + Kk*R*Kk';
    
    % save the diag of the P matrix for covariance stuff 
    errCov(1:3,i) = [P(1); P(5); P(end)];
    
    P_old = P; 
    
end


fig = 1;

% Plot xhat(1,:) with twoSigma1
figure(fig)
twoSigma1 = 2 * sqrt(errCov(1,:));
plot(1:30, xhat(1,:), 'LineWidth', 2); % Bold line for xhat(1,:)
hold on
plot(1:30, twoSigma1, 'r--', 'LineWidth', 1.5); % Dotted line for twoSigma1
plot(1:30, -twoSigma1, 'r--', 'LineWidth', 1.5); % Dotted line for -twoSigma1
hold off
grid on
fig = fig + 1;

% Plot xhat(2,:) with twoSigma2
figure(fig)
twoSigma2 = 2 * sqrt(errCov(2,:));
plot(1:30, xhat(2,:), 'LineWidth', 2); % Bold line for xhat(2,:)
hold on
plot(1:30, twoSigma2, 'r--', 'LineWidth', 1.5); % Dotted line for twoSigma2
plot(1:30, -twoSigma2, 'r--', 'LineWidth', 1.5); % Dotted line for -twoSigma2
hold off
grid on
fig = fig + 1;

% Plot xhat(3,:) with twoSigma3
figure(fig)
twoSigma3 = 2 * sqrt(errCov(3,:));
plot(1:30, xhat(3,:), 'LineWidth', 2); % Bold line for xhat(3,:)
hold on
plot(1:30, twoSigma3, 'r--', 'LineWidth', 1.5); % Dotted line for twoSigma3
plot(1:30, -twoSigma3, 'r--', 'LineWidth', 1.5); % Dotted line for -twoSigma3
hold off
grid on
fig = fig + 1;

% Plot all xhat lines together
figure(fig)
plot(xhat', 'LineWidth', 2); % Bold lines for each xhat
grid on