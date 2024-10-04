%% Homeowrk4 5044

% 10000 random variables
% range a to b

a = -0.5; 
b = 0.5; 
x1 = (b-a).*rand([10000, 1]) + a;
x2 = (b-a).*rand([10000, 1]) + a;

firstPlot = (x1 + x2) / 2;

% PLot with 50 bin histogram 
figure(1)
histogram(firstPlot, 50)
grid on 
xlabel('Random Number')
ylabel('Count')
title('(x1 + x2) / 2')

% now do the second plot 

x3 = (b-a).*rand([10000, 1]) + a;
x4 = (b-a).*rand([10000, 1]) + a;

secondPlot = [x1 + x2 + x3 + x4] / 4; 

figure(2)
histogram(secondPlot, 50)
grid on 
xlabel('Random Number')
ylabel('Count')
title('(x1 + x2 + x3 + x4) / 4')