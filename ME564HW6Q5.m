% ME 564 HW6 Q6

% Define the system dynamics
A = [-3, 2; 0, 1];

% Create a grid for plotting the vector field
[x1, x2] = meshgrid(-5:0.5:5, -5:0.5:5);

% Calculate the vector field
u = A(1, 1) * x1 + A(1, 2) * x2;
v = A(2, 1) * x1 + A(2, 2) * x2;

% Normalize the vectors (for better visualization)
norms = sqrt(u.^2 + v.^2);
u = u ./ norms;
v = v ./ norms;

% Plot the vector field
figure;
quiver(x1, x2, u, v, 0.5);
axis equal;
grid on;
xlabel('x1');
ylabel('x2');
title('Phase Portrait');

% Plot some example trajectories
hold on;
for x0 = [-4:2:4; -4:2:4]
    [t, x] = ode45(@(t, x) A * x, [0, 10], x0);
    plot(x(:, 1), x(:, 2));
end
hold off;