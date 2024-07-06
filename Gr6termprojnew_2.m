% Parameters used for Simulation of Ethanol production
V_i = 1000; % Reactor volume (in liters)
C_glucose_i = 100; % Initial glucose concentration (in g/L)
C_yeast_i = 5; % Initial yeast concentration (in g/L)
C_ethanol_i = 0; % Initial ethanol concentration (in g/L)
mu_max = 0.1; % Maximum specific growth rate (h^-1)
K_s = 25; % Monod constant (g/L)
k_ethanol_production = 0.5; % Ethanol production rate constant (in 1/h)

% Simulation time and step time
Simulation_Time = 6; % Total simulation time (in hours)
Step_Time = 0.05; % Step time (in hours)

% Number of time steps
N = Simulation_Time / Step_Time; % = (Tend - T0)/(step size) =(6 - 0)/(0.05) = 120
 
% Initialization of arrays to store concentrations
C_glucose = zeros(N,1);
C_yeast = zeros(N,1);
C_ethanol = zeros(N,1);
Time = zeros(N,1);

% Setting initial conditions
C_glucose(1,1) = C_glucose_i;
C_yeast(1,1) = C_yeast_i;
C_ethanol(1,1) = C_ethanol_i;
Time(1,1) = 0;

%% Performing simulation
for i = 2:N
    % Calculation of specific growth rate using Monod equation
    mu = mu_max * C_glucose(i-1) / (K_s + C_glucose(i-1)) ;
    % We get the rates of change by calling the function getRatesofChange
    [dC_glucose_dt,dC_yeast_dt,dC_ethanol_dt] = getRatesofChange(C_glucose,C_yeast,mu,k_ethanol_production,i);
    % We get the updated concentrations by calling the explicitEuler function which uses Explicit Euler's method  
    [C_glucose,C_yeast,C_ethanol] = explicitEuler(C_glucose,C_yeast,C_ethanol,dC_glucose_dt,dC_yeast_dt,dC_ethanol_dt,Step_Time,i);
       
    % To ensure concentrations are non-negative
    C_glucose(i) = max(C_glucose(i), 0);
    C_yeast(i) = max(C_yeast(i), 0);
    C_ethanol(i) = max(C_ethanol(i), 0);
    
    % Time updation
    Time(i) = Time(i-1) + Step_Time;
end

% Plotting the graph between concentration of each component and time
figure(1);
plot(Time,C_glucose,'r','LineWidth',1.2);
hold on;
plot(Time,C_yeast,'g','LineWidth',1.2);
plot(Time,C_ethanol,'b','LineWidth',1.2);
xlabel('Time (hours)');
ylabel('Concentration (g/L)');
legend('Glucose','Yeast','Ethanol');
title('Ethanol Production Simulation');
hold off;

% The final concentrations
fprintf('Final concentrations:\n');
fprintf('Glucose: %f g/L\n', C_glucose(end));
fprintf('Yeast: %f g/L\n', C_yeast(end));
fprintf('Ethanol: %f g/L\n', C_ethanol(end));

% Function to calculate the rate of change of each component
function [dC_glucose_dt,dC_yeast_dt,dC_ethanol_dt] = getRatesofChange(C_glucose,C_yeast,mu,k_ethanol_production,j)

       dC_glucose_dt = -mu * C_glucose(j-1) * C_yeast(j-1);
       dC_yeast_dt = mu * C_glucose(j-1) * C_yeast(j-1) - k_ethanol_production * C_yeast(j-1);
       dC_ethanol_dt = k_ethanol_production * mu * C_glucose(j-1) * C_yeast(j-1);
       
end

% Updating concentrations by using Explicit Euler's method 
% Y(i+1) = Y(i) + h*f( Xi , Yi )
function [C_glucose,C_yeast,C_ethanol] = explicitEuler(C_glucose,C_yeast,C_ethanol,f_glucose,f_yeast,f_ethanol,h,k)

     C_glucose(k) = C_glucose(k-1) + h * f_glucose ;
     C_yeast(k) = C_yeast(k-1) + h * f_yeast ; 
     C_ethanol(k) = C_ethanol(k-1) + h * f_ethanol ;

end
