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

% stm multiplied by ICs - U is zero so no worry about it!

x_k = F * x0; 

for k = 1:length(t)-1
    
    x_k_1{k+1} = F * x_k;
    
    x_k = x_k_1{k+1}; 
end

x_k_1{1} = F * x0;

Stateplot = cell2mat(x_k_1)';

% Subplot 1
subplot(4, 1, 1)
plot(10*(1:length(Stateplot(:, 1))), Stateplot(:, 1), 'LineWidth', 1.5) % Increase line width
grid on
ylabel('$\Delta r$', 'Interpreter', 'latex', 'FontSize', 12)


% Subplot 2
subplot(4, 1, 2)
plot(10*(1:length(Stateplot(:, 1))), Stateplot(:, 2), 'LineWidth', 1.5) % Increase line width
grid on
ylabel('$\Delta \dot{r}$', 'Interpreter', 'latex', 'FontSize', 12)


% Subplot 3
subplot(4, 1, 3)
plot(10*(1:length(Stateplot(:, 1))), Stateplot(:, 3), 'LineWidth', 1.5) % Increase line width
ylabel('$\Delta \theta$', 'Interpreter', 'latex', 'FontSize', 12)
grid on


% Subplot 4
subplot(4, 1, 4)
plot(10*(1:length(Stateplot(:, 1))), Stateplot(:, 4), 'LineWidth', 1.5) % Increase line width
grid on
ylabel('$\Delta \dot{\theta}$', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('Time (sec)', 'FontSize', 12) % Increase x-label font size


sgtitle('Propagation of Perturbation States')


% Now total state

% We have been given the nominal state: 

for i = 1:length(t)
    state_nominal{i} = [r0; 0; w0*t(i)+sqrt(398600/(r0^3)); w0];
end

% overwite first time step with appriopirate values
state_nominal{1} = [r0; 0; 0; w0];

nom_state = cell2mat(state_nominal)';

total_state = nom_state + Stateplot;

figure(101)

subplot(4,1,1)
plot(10*(1:length(total_state(:,1))), total_state(:,1), 'LineWidth', 1.5)
grid on
ylabel('$r$', 'Interpreter', 'latex', 'FontSize', 12)

subplot(4,1,2)
plot(10*(1:length(total_state(:,1))), total_state(:,2), 'LineWidth', 1.5)
grid on
ylabel('$\dot{r}$', 'Interpreter', 'latex', 'FontSize', 12)

subplot(4,1,3)
plot(10*(1:length(total_state(:,1))), total_state(:,3), 'LineWidth', 1.5)
grid on
ylabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 12)

subplot(4,1,4)
plot(10*(1:length(total_state(:,1))), total_state(:,4), 'LineWidth', 1.5)
grid on
ylabel('$\dot{\theta}$', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('Time (sec)', 'FontSize', 12)

sgtitle('Total State Over Orbit')


%% part b - ODE45 
% ODE45 

initCond = x0 + [r0; 0; 0; w0];
tspan = 0:10:5400;
options = odeset('RelTol', 1e-5);
[t,Y] = ode45(@myOde, tspan, initCond, options);

figure(4)
set(gcf, 'Position', [100, 100, 800, 600]); % Resize figure

% Plot 1
subplot(4,1,1)
plot(t, Y(:,1), 'LineWidth', 2);
grid on
ylabel('r', 'FontSize', 12)
title('Dynamics of the System', 'FontSize', 14)
xlim([min(t) max(t)])
set(gca, 'FontSize', 12)

% Plot 2
subplot(4,1,2)
plot(t, Y(:,2), 'LineWidth', 2);
grid on
ylabel('$\dot{r}$', 'Interpreter', 'latex', 'FontSize', 12);
xlim([min(t) max(t)])
set(gca, 'FontSize', 12)

% Plot 3
subplot(4,1,3)
plot(t, Y(:,3), 'LineWidth', 2);
grid on
ylabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 12);
xlim([min(t) max(t)])
set(gca, 'FontSize', 12)

% Plot 4
subplot(4,1,4)
plot(t, Y(:,4), 'LineWidth', 2);
grid on
ylabel('$\dot{\theta}$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12) % Add x-label for the last subplot
xlim([min(t) max(t)])
set(gca, 'FontSize', 12)

% Adjust the layout
sgtitle('State Variables Over Time ODE45 Integration', 'FontSize', 16) % Super title for the entire figure


% Difference between nominal and actual trajectories 

diff_nom_act = Y - nom_state; 

figure(102)
set(gcf, 'Position', [100, 100, 800, 600]); % Resize figure

subplot(4,1,1)
plot(10*(1:length(diff_nom_act(:,1))), diff_nom_act(:,1), 'LineWidth', 2)
grid on
ylabel('$\dot{r}$', 'Interpreter', 'latex', 'FontSize', 12);

subplot(4,1,2)
plot(10*(1:length(diff_nom_act(:,1))), diff_nom_act(:,2), 'LineWidth', 2)
grid on
ylabel('$\dot{r}$', 'Interpreter', 'latex', 'FontSize', 12);

subplot(4,1,3)
plot(10*(1:length(diff_nom_act(:,1))), diff_nom_act(:,3), 'LineWidth', 2)
grid on
ylabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 12);

subplot(4,1,4)
plot(10*(1:length(diff_nom_act(:,1))), diff_nom_act(:,4), 'LineWidth', 2)
grid on
ylabel('$\dot{\theta}$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time (sec)')

sgtitle('Perturbation States using ODE45', 'FontSize', 16) % Super title for the entire figure



% compare ODE 45 integration to discretized integration

% discretized perturbations
% Stateplot

% ode45 perutbations
% diff_nom_act

dt_ode_diff = diff_nom_act - Stateplot; 

figure(105)
set(gcf, 'Position', [100, 100, 800, 600]); % Resize figure

subplot(4,1,1)
plot(10*(1:length(dt_ode_diff(:,1))), dt_ode_diff(:,1), 'LineWidth', 2)
grid on
ylabel('$\dot{r}$', 'Interpreter', 'latex', 'FontSize', 12);

subplot(4,1,2)
plot(10*(1:length(dt_ode_diff(:,1))), dt_ode_diff(:,2), 'LineWidth', 2)
grid on
ylabel('$\dot{r}$', 'Interpreter', 'latex', 'FontSize', 12);

subplot(4,1,3)
plot(10*(1:length(dt_ode_diff(:,1))), dt_ode_diff(:,3), 'LineWidth', 2)
grid on
ylabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 12);

subplot(4,1,4)
plot(10*(1:length(dt_ode_diff(:,1))), dt_ode_diff(:,4), 'LineWidth', 2)
grid on
ylabel('$\dot{\theta}$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time (s)')

sgtitle('Difference in DT and ODE perturbation states', 'FontSize', 16) % Super title for the entire figure


% difference in total state for DT and ode45

% DT total state
% total_state

% ode45 total state
%Y


dt_ode_diff_total = Y - total_state;

figure(107)
set(gcf, 'Position', [100, 100, 800, 600]); % Resize figure

subplot(4,1,1)
plot(10*(1:length(dt_ode_diff_total(:,1))), dt_ode_diff_total(:,1), 'LineWidth', 2)
grid on
ylabel('$\dot{r}$', 'Interpreter', 'latex', 'FontSize', 12);

subplot(4,1,2)
plot(10*(1:length(dt_ode_diff_total(:,1))), dt_ode_diff_total(:,2), 'LineWidth', 2)
grid on
ylabel('$\dot{r}$', 'Interpreter', 'latex', 'FontSize', 12);

subplot(4,1,3)
plot(10*(1:length(dt_ode_diff_total(:,1))), dt_ode_diff_total(:,3), 'LineWidth', 2)
grid on
ylabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 12);

subplot(4,1,4)
plot(10*(1:length(dt_ode_diff_total(:,1))), dt_ode_diff_total(:,4), 'LineWidth', 2)
grid on
ylabel('$\dot{\theta}$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time (s)')

sgtitle('Difference in DT and ODE total states', 'FontSize', 16) % Super title for the entire figure


%% Symbolic fun

syms s
syms t

to = [0 1 0; 0 0 1; 0 0 0] - eye(3)*s

det(to)


[1 -1 1; 0 0 0; 0 0 0]*noel*inv([1 -1 1; 0 0 0; 0 0 0])