% ME 564 HW6 Q3

% Given matrices A and B
A = [-2.6, 0.25, -38, 0;
     -0.075, -0.27, 4.4, 0;
     0.078, -0.99, -0.23, 0.052;
     1.0, 0.078, 0, 0];
     
B = [17, 7;
     0.82, -3.2;
     0, 0.046;
     0, 0];

% a
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

% Create state-space model for simulation
sys = ss(A, B, eye(4), zeros(4, 2));

% Time vector
t = linspace(0, 5, 500);

% b
% Input signals for u1 active
u1 = zeros(size(t));
u1(t >= 0 & t <= 1) = 1;
u2 = zeros(size(t));
u = [u1; u2]'; % Transpose to make it 500x2
% u = [u1';u2'];

% Simulate the system
[y, t] = lsim(sys, u, t);

% Plot the states
figure;
hold on;
title('Dynamical Response for u1 active');
plotStates(t, y);

% c
% Input signals for u2 active
u2 = zeros(size(t));
u2(t >= 0 & t <= 1) = 1;
u1 = zeros(size(t));
u = [u1'; u2'];  % Transpose to make it 500x2

% Simulate the system
[y, t] = lsim(sys, u, t);

% Plot the states
figure;
hold on;
title('Dynamical Response for u2 active');
plotStates(t, y);

% Function to plot states
function plotStates(t, y)
    plot(t, y(:, 1), 'r', 'DisplayName', 'x1: Roll rate');
    plot(t, y(:, 2), 'g', 'DisplayName', 'x2: Yaw rate');
    plot(t, y(:, 3), 'b', 'DisplayName', 'x3: Sideslip angle');
    plot(t, y(:, 4), 'm', 'DisplayName', 'x4: Roll attitude');
    xlabel('Time (s)');
    ylabel('State values');
    legend;
    hold off;
end
