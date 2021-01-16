%% Define parameters
transmission_rate = 1.4; % Transmission rate 
clinical_infection = 0.4; % Probability that the infection becomes clinical
incubation_period = 4.5; % Average incubation period
infection_period = 5; % Average infection period
average_annual_salary = 580000; % Average annual salary in NOK
k = 1:5; % Defining the different clusters
%C = zeros(size(k, 2), size(k, 2)); % Initialize contact matrix C
C = [100 20 30 10 10; 30 100 30 5 5; 20 30 100 10 10; 5 10 30 100 20; 5 5 20 30 100]; % Contact matrix
%% Defining functions
infectious_probability = 1 - exp(-1/incubation_period);
recovery_probability = 1 - exp(-1/infection_period);
%% Initial Conditions S1 I1 E1 R1 S2 I2  
IC = [100 0 0 0 100 0 0 0 100 0 0 0 100 0 0 0 100 0 0 0];
%% Time span
tinitial = 0; % Starting at day = 0
tfinal = 500; % Modelling for the next 500 days
timespan = tinitial:tfinal;
%% Model/ system of differential equations
% Define functions
syms S1(t) S2(t) S3(t) S4(t) S5(t)
syms I1(t) I2(t) I3(t) I4(t) I5(t)
syms E1(t) E2(t) E3(t) E4(t) E5(t)
syms R1(t) R2(t) R3(t) R4(t) R5(t)
% Define ODEs
ode1 = diff(S1,t) == S1 - transmission_rate.*S1.*(C(1,1).*I1+C(1,2).*I2+C(1,3).*I3+C(1,4).*I4+C(1,5).*I5);
ode2 = diff(E1,t) == (1 - infectious_probability).*E1 + transmission_rate.*S1.*(C(1,1).*I1+C(1,2).*I2+C(1,3).*I3+C(1,4).*I4+C(1,5).*I5);
ode3 = diff(I1,t) == clinical_infection.*infectious_probability.*E1 + (1 - recovery_probability).*I1;
ode4 = diff(R1,t) == R1 + recovery_probability.*I1;
ode5 = diff(S2,t) == S2 - transmission_rate.*S2.*(C(2,1).*I1+C(2,2).*I2+C(2,3).*I3+C(2,4).*I4+C(2,5).*I5);
ode6 = diff(E2,t) == (1 - infectious_probability).*E2 + transmission_rate.*S2.*(C(2,1).*I1+C(2,2).*I2+C(2,3).*I3+C(2,4).*I4+C(2,5).*I5);
ode7 = diff(I2,t) == clinical_infection.*infectious_probability.*E2 + (1 - recovery_probability).*I2;
ode8 = diff(R2,t) == R2 + recovery_probability.*I2;
ode9 = diff(S3,t) == S3 - transmission_rate.*S3.*(C(3,1).*I1+C(3,2).*I2+C(3,3).*I3+C(3,4).*I4+C(3,5).*I5);
ode10 = diff(E3,t) == (1 - infectious_probability).*E3 + transmission_rate.*S3.*(C(3,1).*I1+C(3,2).*I2+C(3,3).*I3+C(3,4).*I4+C(3,5).*I5);
ode11 = diff(I3,t) == clinical_infection.*infectious_probability.*E3 + (1 - recovery_probability).*I3;
ode12 = diff(R3,t) == R3 + recovery_probability.*I3;
ode13 = diff(S4,t) == S4 - transmission_rate.*S4.*(C(4,1).*I1+C(4,2).*I2+C(4,3).*I3+C(4,4).*I4+C(4,5).*I5);
ode14 = diff(E4,t) == (1 - infectious_probability).*E4 + transmission_rate.*S4.*(C(4,1).*I1+C(4,2).*I2+C(4,3).*I3+C(4,4).*I4+C(4,5).*I5);
ode15 = diff(I4,t) == clinical_infection.*infectious_probability.*E4 + (1 - recovery_probability).*I4;
ode16 = diff(R4,t) == R4 + recovery_probability.*I4;
ode17 = diff(S5,t) == S5 - transmission_rate.*S5.*(C(5,1).*I1+C(5,2).*I2+C(5,3).*I3+C(5,4).*I4+C(5,5).*I5);
ode18 = diff(E5,t) == (1 - infectious_probability).*E5 + transmission_rate.*S5.*(C(5,1).*I1+C(5,2).*I2+C(5,3).*I3+C(5,4).*I4+C(5,5).*I5);
ode19 = diff(I5,t) == clinical_infection.*infectious_probability.*E5 + (1 - recovery_probability).*I5;
ode20 = diff(R5,t) == R5 + recovery_probability.*I5;
odes = [ode1; ode2; ode3; ode4; ode5; ode6; ode7; ode8; ode9; ode10; ode11; ode12; ode13; ode14; ode15; ode16; ode17; ode18; ode19; ode20];
[V, Y] = odeToVectorField(odes);
M = matlabFunction(V, 'vars', {'t', 'Y'});
%% Solution
sol = ode23s(M, timespan, IC); 