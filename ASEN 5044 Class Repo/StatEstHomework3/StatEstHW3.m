%% 5044 Homework 3

% Problem 1

A = [0 1 0 0; -2 0 1 0; 0 0 0 1; 1 0 -2 0]; 
B = [0 0; -1 0; 0 0; 1 1]; 
C = [1 0 0 0; 0 1 0 -1]; 
D = [0 0; 0 0];

deltaT = 0.05; 

% part a 
Ahat = [[A], [B]; zeros(2,6)]; 

Ahatexp = expm(Ahat*deltaT); 

F = Ahatexp(1:4, 1:4)
G = Ahatexp(1:4, 5:6)

H = C; 
M = D; 

% smapling rate compare to the Nyquist rate? 
[eVectors, eValues] = eig(A); 

sampleRate = 2*pi / deltaT

NyquitCriteria = pi / norm(eValues(1))

if deltaT < NyquitCriteria
    disp('Meets the Nyquist threshold!')
else
    disp('DOES NOT MEET NYQUIST')
end

% part b

observabilityMatrix = [H; H*F; H*F^2; H*F^3]; 

if rank(observabilityMatrix) == 4
    disp('Entire system observable!')
else
    disp('Not observable')
end

% part c
load('hw3problem1data')

[x0] = InputOutputObservability(Udata, Ydata, F, G, H, M)

% tranpose for later on ease!
Udata = Udata'; 
Ydata = Ydata';

% part d

t = 0:deltaT:5; 

x_k = x0; 

for k = 1:length(t)
    % states over time 
    x_future{k+1} = F*x_k + G*Udata(:,k);
    
    if k == 1
        % do nothing on first time because time = 0 
    else
        % measured output over time
        measuredOutput{k-1} = H*x_k + M*Udata(:,k);
    end
    
    % for plotting the state
    CurentStateHist{k} = x_k; 
    
    % update state for next go around
    x_k = x_future{k+1};
end

CurentStateHist = cell2mat(CurentStateHist);
CurentStateHist = CurentStateHist'; 

% plot the states over time
fig = 1;

figure(fig)
plot(t, CurentStateHist)
grid on
xlabel('Time [sec]')
ylabel('System States')
legend('$x_1$', '$x_2$', '$x_3$', '$x_4$', 'Interpreter', 'latex')
fig = fig + 1; 
title('All System States Over Time')


% plot the measured output over time 
measuredOutput = cell2mat(measuredOutput);
measuredOutput = measuredOutput';

figure(fig)
plot(t(2):0.05:t(end), measuredOutput)
grid on
xlabel('Time [sec]')
title('Recorded Measured Output')
legend('$y_1$', '$y_2$', 'Interpreter', 'latex')
fig = fig + 1; 


figure(fig)
plot(t(2):deltaT:t(end), Ydata)
grid on
xlabel('Time [sec]')
title('Given Ydata')
legend('$y_1$', '$y_2$', 'Interpreter', 'latex')
fig = fig + 1; 


% plot the differences in the recorded and predicted

y1_diff = Ydata(1,:) - measuredOutput(:,1)';
y2_diff = Ydata(2,:) - measuredOutput(:,2)';

figure(fig)
plot(t(2):deltaT:t(end), y1_diff, t(2):deltaT:t(end), y2_diff)
grid on 
title('diff')
fig = fig+1; 


%% Problem 2 

