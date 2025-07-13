clear all
close all

% Define Parameters
m1 = 1;
m2 = 1;
I1 = 1;
I2 = 1;
l1 = 1;
g = 0;

% Time span
tspan = [0 20];

% Initial condition (equilibrium point)
x0 = [0; 0; 0; 0];

% Solve nonlinear system
[tn, xn] = ode45(@(t, x) nonlinear_system(t, x, m1, m2, I1, I2, l1, g), tspan, x0);

% Solve linear system
[t, xl] = ode45(@(t, x) linear_system(t, x, m1, m2, I1, I2, l1, g), tspan, x0);

% Plot Results
subplot(2,1,1);
plot(tn, xn(:,1), '-', t, xl(:,1), '-.');
title('Theta1 vs Time');
xlabel('Time (s)');
ylabel('Theta1 (rad)');

subplot(2,1,2);
plot(tn, xn(:,3), '-', t, xl(:,3), '-.');
title('d2 vs Time');
xlabel('Time (s)');
ylabel('d2 (m)');

function dx = nonlinear_system(t, x, m1, m2, I1, I2, l1, g)
    theta1 = x(1);
    theta1_dot = x(2);
    d2 = x(3);
    d2_dot = x(4);
    
    tau1 = 0.1 * sin(t);
    tau2 = 0.1 * sin(t);
    
    theta1_ddot = (tau1 - 2*m2*d2*theta1_dot*d2_dot - (m1*l1 + m2*d2)*g*cos(theta1)) / (m1*l1^2 + I1 + I2 + m2*d2^2);
    d2_ddot = (tau2 + m2*d2*theta1_dot^2 - m2*g*sin(theta1)) / m2;
    
    dx = [theta1_dot; theta1_ddot; d2_dot; d2_ddot];
end

function dx_tilda = linear_system(t, x_tilda, m1, m2, I1, I2, l1, g)
    % Define A, B, and c matrices
    A = [0, 1, 0, 0; 
         0, -g, 0, 0; 
         0, 0, 0, 1; 
         0, 0, 0, 0];
    
    B = [0, 0; 
         1/(m1*l1^2 + I1 + I2 + g*m2), 0; 
         0, 0; 
         0, 1/m2];
    
    c = [(m1*l1 + 3*m2)*g; 0; 0; 0];
    
    % Define u_tilda (deviation from the operating point)
    % Assuming u_op = 0 for simplicity
    u_tilda = [0.1 * sin(t); 0.1 * sin(t)] - [0; 0];
   
    % Compute dx_tilda
    dx_tilda = A * x_tilda + B * u_tilda + c;
end
