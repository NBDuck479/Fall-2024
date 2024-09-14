%% HW 1

% Problem 3

Iy = 750; 
Iz = 1000; 
Ix = 500; 
Po = 20; 
deltaT = 0.1; 

A = [0 0 0; 0 0 Po*(Ix - Iz)/Iy; 0 Po*(Iy - Ix)/Iz 0]

expm(A * deltaT);
%% 
deltaq0 = 0.1; 
deltap0 = 0; 
deltar0 = 0; 
timeSeries = 0:deltaT:5; 

for i = 1:length(timeSeries)
    State{i} = expm(A * timeSeries(i)) * [deltap0; deltaq0; deltar0];
end

StateVector = cell2mat(State)';

figure(1)
subplot(3,1,1)
plot(timeSeries, StateVector(:,1))
xlabel('Time (sec)')
ylabel('\Delta p')
grid on

subplot(3,1,2)
plot(timeSeries, StateVector(:,2))
xlabel('Time (sec)')
ylabel('\Delta q')
grid on

subplot(3,1,3)
plot(timeSeries, StateVector(:,3))
xlabel('Time (sec)')
ylabel('\Delta r')
grid on