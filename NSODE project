clc;
clear;
close all;

%% -------------------- Parameters --------------------
A  = 1449401;     % Total population
r  = 0.5;         % Transmission rate (per month)
a  = 0.111111;    % Recovery rate (per month)
mu = 0.001167;    % Death rate (per month)

%% -------------------- Time Settings --------------------
t0 = 0;           % Initial time (months)
tf = 50;          % Final time (months)
h  = 0.1;         % Step size
t  = t0:h:tf;
n  = length(t);

%% -------------------- Initial Conditions --------------------
S = zeros(1,n);
I = zeros(1,n);
R = zeros(1,n);

S(1) = 1446093;    % Susceptible
I(1) = 1885;       % Infected
R(1) = 1423;       % Recovered

%% -------------------- RK4 Method --------------------
for k = 1:n-1
    
    N = S(k) + I(k) + R(k);
    
    % k1
    k1S = A - (r*S(k)*I(k))/N - mu*S(k);
    k1I = (r*S(k)*I(k))/N - (mu + a)*I(k);
    k1R = a*I(k) - mu*R(k);
    
    % k2
    S2 = S(k) + h*k1S/2;
    I2 = I(k) + h*k1I/2;
    R2 = R(k) + h*k1R/2;
    N2 = S2 + I2 + R2;
    
    k2S = A - (r*S2*I2)/N2 - mu*S2;
    k2I = (r*S2*I2)/N2 - (mu + a)*I2;
    k2R = a*I2 - mu*R2;
    
    % k3
    S3 = S(k) + h*k2S/2;
    I3 = I(k) + h*k2I/2;
    R3 = R(k) + h*k2R/2;
    N3 = S3 + I3 + R3;
    
    k3S = A - (r*S3*I3)/N3 - mu*S3;
    k3I = (r*S3*I3)/N3 - (mu + a)*I3;
    k3R = a*I3 - mu*R3;
    
    % k4
    S4 = S(k) + h*k3S;
    I4 = I(k) + h*k3I;
    R4 = R(k) + h*k3R;
    N4 = S4 + I4 + R4;
    
    k4S = A - (r*S4*I4)/N4 - mu*S4;
    k4I = (r*S4*I4)/N4 - (mu + a)*I4;
    k4R = a*I4 - mu*R4;
    
    % Update
    S(k+1) = S(k) + (h/6)*(k1S + 2*k2S + 2*k3S + k4S);
    I(k+1) = I(k) + (h/6)*(k1I + 2*k2I + 2*k3I + k4I);
    R(k+1) = R(k) + (h/6)*(k1R + 2*k2R + 2*k3R + k4R);
end

%% -------------------- Plot Results --------------------
figure;
plot(t,S,'b','LineWidth',2); hold on;
plot(t,I,'r','LineWidth',2);
plot(t,R,'g','LineWidth',2);
xlabel('Time (months)');
ylabel('Population');
legend('Susceptible','Infected','Recovered');
title('SIR Model for Tuberculosis using RK4 Method');
grid on;
