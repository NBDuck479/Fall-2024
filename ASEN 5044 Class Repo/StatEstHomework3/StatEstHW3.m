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

F = Ahatexp(1:4, 1:4);
G = Ahatexp(1:4, 4:5); 

H = C; 
M = D; 

% smapling rate compare to the Nyquist rate? 
[eVectors, eValues] = eig(A); 

sampleRate = 2*pi / deltaT; 

NyquitCriteria = pi / norm(eValues(1));

if deltaT < NyquitCriteria
    disp('Meets the Nyquist threshold!')
else
    disp('DOES NOT MEET NYQUIST')
end