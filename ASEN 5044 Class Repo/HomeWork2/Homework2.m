% Homework 2

% Part c
k = 398600; 
r0 = 6678; 
deltaT = 10; 
w0 = sqrt(k/(r0^3))

% plug in give values for A 

A = [0 1 0 0; w0^2+2*k/(r0^3) 0 0 2*r0*w0; 0 0 0 1; 0 -2*w0/r0 0 0]

B = [0 0; 1 0; 0 0; 0 1/r0]

C = [1 0 0 0; 0 0 1 0]
D = [0 0; 0 0]

% convert to discrete time
F = expm(A*deltaT)

% augemented state
Ahat = [[A],[B]; zeros(2,6)]
augmented = expm(Ahat*deltaT)



%% problem 4
% initial conditions r rdot theta thetadot
x0 = [10; -.5; 0; 2.5*10^-5];
t = 0:10:5400;
u = zeros(2, length(t));

% create sys for lsim plotting of discrete time 
G = augmented(1:4, 5:6)
H = C
M = D

sys = ss(F,G,H,D)

lsim(sys,u,t,x0)

% stm multiplied by ICs - U is zero so no worry about it!

x_k = F * x0; 

for k = 1:length(t)
    
    x_k_1{k+1} = F * x_k;
    
    x_k = x_k_1{k+1}; 
end

x_k_1{1} = F * x0;

Stateplot = cell2mat(x_k_1)'

% Subplot 1
subplot(4, 1, 1)
plot(Stateplot(:, 1), 'LineWidth', 1.5) % Increase line width
grid on
ylabel('$\Delta r$', 'Interpreter', 'latex', 'FontSize', 12)
title('Change in Position', 'FontSize', 14) % Add title
xlim([1, length(Stateplot(:, 1))]) % Set x-axis limits

% Subplot 2
subplot(4, 1, 2)
plot(Stateplot(:, 2), 'LineWidth', 1.5) % Increase line width
grid on
ylabel('$\Delta \dot{r}$', 'Interpreter', 'latex', 'FontSize', 12)
title('Change in Velocity', 'FontSize', 14) % Add title
xlim([1, length(Stateplot(:, 2))]) % Set x-axis limits

% Subplot 3
subplot(4, 1, 3)
plot(Stateplot(:, 3), 'LineWidth', 1.5) % Increase line width
ylabel('$\Delta \theta$', 'Interpreter', 'latex', 'FontSize', 12)
grid on
title('Change in Angle', 'FontSize', 14) % Add title
xlim([1, length(Stateplot(:, 3))]) % Set x-axis limits

% Subplot 4
subplot(4, 1, 4)
plot(Stateplot(:, 4), 'LineWidth', 1.5) % Increase line width
grid on
ylabel('$\Delta \dot{\theta}$', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('Time (sec)', 'FontSize', 12) % Increase x-label font size
title('Change in Angular Velocity', 'FontSize', 14) % Add title
xlim([1, length(Stateplot(:, 4))]) % Set x-axis limits


% Now get the total state 

x_k = F * x0; 

for k = 1:length(t)
    
    x_k_1{k+1} = F * x_k;
    
    x_k = x_k_1{k+1}; 
end

x_k_1{1} = F * x0;



% part b 
% ODE45 
tspan = [0 100]; 
[t,Y] = ode45(@myOde, tspan, x0)

figure(2)
subplot(4,1,1)
plot(Y(:,1))
grid on
ylabel('\Delta r')

subplot(4,1,2)
plot(Y(:,2))
grid on
ylabel('$\Delta \dot{r}$', 'Interpreter', 'latex');

subplot(4,1,3)
plot(Y(:,3))
grid on
ylabel('$\Delta \theta$', 'Interpreter', 'latex');

subplot(4,1,4)
plot(Y(:,4))
grid on
ylabel('$\Delta \dot{\theta}$', 'Interpreter', 'latex');