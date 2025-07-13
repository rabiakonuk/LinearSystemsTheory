% ME 564 HW6 Q4

% a
% Given values
M = 0.5;  % Mass of the cart (kg)
m = 0.2;  % Mass of the pendulum (kg)
b = 0.1;  % Friction coefficient
I = 0.006;  % Moment of inertia (kg*m^2)
g = 9.8;  % Acceleration due to gravity (m/s^2)
l = 0.3;  % Length of the pendulum (m)

% Calculate the denominator for the A and B matrices
denominator = I * (M + m) + M * m * l^2;

% Define the A matrix
A = [0, 1, 0, 0;
     0, -((I + m * l^2) * b) / denominator, (m^2 * g * l^2) / denominator, 0;
     0, 0, 0, 1;
     0, -(m * l * b) / denominator, (m * g * l * (M + m)) / denominator, 0];

% Define the B matrix
B = [0;
     (I + m * l^2) / denominator;
     0;
     (m * l) / denominator];

% Display the A and B matrices
disp('Matrix A:');
disp(A);
disp('Matrix B:');
disp(B);

% b
% Calculate the controllability matrix C
n = size(A, 1);  % Number of states
C = B;
for i = 1:n-1
    C = [C, A^i * B];
end

% Check if the system is controllable
rank_C = rank(C);
if rank_C == n
    disp('The system is controllable.');
else
    disp('The system is not controllable.');
end

% c
% Given time t
t = 0.1;

% Compute the state-transition matrix Phi(t, 0)
Phi_t_0 = expm(A * t);

% Display the state-transition matrix
disp('State-transition matrix Phi(t, 0) for t = 0.1:');
disp(Phi_t_0);

% d
% Given time t
t = 1;

% Compute the state-transition matrix Phi(t, 0)
Phi_t_0 = expm(A * t);

% Display the state-transition matrix
disp('State-transition matrix Phi(t, 0) for t = 1:');
disp(Phi_t_0);

% e
% Create state-space model
sys = ss(A, B, eye(4), zeros(4, 1));

% Time vector
t = linspace(0, 10, 1000);

% Unit step input F(t) = 1 for t >= 0
F = ones(size(t));

% Initial state z(0) = 0
z0 = zeros(4, 1);

% Simulate the system
[y, t] = lsim(sys, F, t, z0);

% Plot x and phi on the same graph
figure;
hold on;
plot(t, y(:, 1), 'r', 'DisplayName', 'x: Deviation from equilibrium');
plot(t, y(:, 3), 'b', 'DisplayName', '\phi: Deviation from vertical');
xlabel('Time (s)');
ylabel('State values');
title('Response to Unit Step Input F(t) = 1');
legend;
hold off;
