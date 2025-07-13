% ME 564 HW 9 
% q3

% a
% 1: Find Eigenvalues and Eigenvectors
A = [-2 -1; 0 -3];
[V, D] = eig(A);  % V: eigenvectors, D: eigenvalues as a diagonal matrix

disp('Eigenvalues:');
disp(diag(D));

disp('Eigenvectors:');
disp(V);

% 2: Form the Transformation Matrix V
% V is already obtained from eig function

% 3: Transform the System to Modal Form
A_modal = inv(V) * A * V;
B_modal = inv(V) * [2; 1];

disp('Modal Form Matrices:');
disp('A_modal:');
disp(A_modal);

disp('B_modal:');
disp(B_modal);


% b 
% Given matrices
A_modal = [-2 0; 0 -3];
B_modal = [2; 1];

% Time vector
% assume:
t = linspace(0, 5, 100); 

% Define the input function u(t)
% assume input function is:
u = @(t) sin(t);  

% Compute Lu(t) - Zero-state response
Lu = zeros(size(B_modal, 1), length(t));
for i = 1:length(t)
    Lu(:, i) = integral(@(tau) expm(A_modal*(t(i)-tau))*B_modal*u(tau), 0, t(i), 'ArrayValued', true);
end

% Compute Lo(t) - Zero-input response with x(0) = 0
Lo = zeros(size(A_modal, 1), length(t));
for i = 1:length(t)
    Lo(:, i) = expm(A_modal*t(i)) * zeros(size(A_modal, 1), 1);
end

% Display the results
figure;

subplot(2, 1, 1);
plot(t, Lu(1, :), 'r', t, Lu(2, :), 'b', 'LineWidth', 2);
title('Zero-State Response (Lu)');
xlabel('Time');
ylabel('State Variables');
legend('Lu (z1)', 'Lu (z2)');

subplot(2, 1, 2);
plot(t, Lo(1, :), 'r', t, Lo(2, :), 'b', 'LineWidth', 2);
title('Zero-Input Response (Lo)');
xlabel('Time');
ylabel('State Variables');
legend('Lo (z1)', 'Lo (z2)');


% c

% Given matrices
A = [-2 -1; 0 -3];
B = [2 1];
C = [1; -1]; % Make C a column vector for proper transpose

% Step 1: Find Adjoint Matrix A*
A_star = conj(A.');

u = @(t) sin(t);


% Step 2: Compute Pu(t) - Zero-state response for adjoint system
Pu = zeros(size(A_star, 1), length(t));
for i = 1:length(t)
    % Define the matrix exponential function
    expAt = @(t) expm(A_star * t);

    % Define the system dynamics for the integral
    % dynamics = @(tau, z) expAt(tau) * C.' * u(tau);
    
    % expAt_t = expAt(t(1));  
    % C_t = C.' * u(t(1)).';
    % disp(size(expAt_t));
    % disp(size(C_t));

    dynamics = @(tau, z) expAt(tau) * (C.' * u(tau)).';
    
    % Initialize state variable for Pu
    z0 = zeros(size(A_star, 2), 1);

    % Integrate the system using a loop
    for j = 1:length(t)
        [~, z] = ode45(@(tau, z) dynamics(tau, z), [eps t(j)], z0); % Start from a small positive value
        z0 = z(end, :).';  % Update initial condition for the next step
    end

    % Store the final state for Pu
    Pu(:, i) = z0;
end

% Step 3: Compute Po(t) - Zero-input response for adjoint system
Po = zeros(size(A_star, 1), length(t));
x0 = [0; 0]; % Initial conditions for the adjoint system
for i = 1:length(t)
    % Initialize state variable for Po
    z0 = x0;

    % Integrate the system using a loop
    for j = 1:length(t)
        [~, z] = ode45(@(tau, z) A_star * z, [eps t(j)], z0); % Start from a small positive value
        z0 = z(end, :).';  % Update initial condition for the next step
    end

    % Store the final state for Po
    Po(:, i) = z0;
end

% Display the results
figure;

subplot(2, 1, 1);
plot(t, Pu, 'k', 'LineWidth', 2);
title('Zero-State Response (Pu) - Adjoint System');
xlabel('Time');
ylabel('State Variable');
legend('Pu');

subplot(2, 1, 2);
plot(t, Po, 'k', 'LineWidth', 2);
title('Zero-Input Response (Po) - Adjoint System');
xlabel('Time');
ylabel('State Variable');
legend('Po');

% d

% Calculate the adjoint matrices for Lo_star
A_star_lo = conj(A_modal.');  % Adjoint of A_modal
B_star_lo = V * B_modal;       % Adjoint of B_modal
C_star_lo = C.' * inv(V);      % Adjoint of C

% Time vector for response comparison
t_compare = linspace(0, 5, 100);

% Compute Lo_star(t) - Zero-input response for adjoint system
Lo_star = zeros(size(A_star_lo, 1), length(t_compare));
x0_lo_star = [0; 0]; % Initial conditions for Lo_star
for i = 1:length(t_compare)
    % Initialize state variable for Lo_star
    z0_lo_star = x0_lo_star;

    % Integrate the system using a loop
    for j = 1:length(t_compare)
        [~, z_lo_star] = ode45(@(tau, z) A_star_lo * z, [eps t_compare(j)], z0_lo_star); % Start from a small positive value
        z0_lo_star = z_lo_star(end, :).';  % Update initial condition for the next step
    end

    % Store the final state for Lo_star
    Lo_star(:, i) = z0_lo_star;
end

% Compare Lo_star with Pu
figure;

subplot(2, 1, 1);
plot(t_compare, Pu, 'k', 'LineWidth', 2);
title('Zero-State Response (Pu)');
xlabel('Time');
ylabel('State Variable');
legend('Pu');

subplot(2, 1, 2);
plot(t_compare, Lo_star, 'k--', 'LineWidth', 2);
title('Zero-Input Response (Lo_star)');
xlabel('Time');
ylabel('State Variable');
legend('Lo\_star');

% Display legend for comparison
legend({'Pu', 'Lo\_star'});


% % e 
% % Given matrices
% A = [-2 -1; 0 -3];
% B = [2 1];
% C = [1; -1]; % Make C a column vector for proper transpose
% 
% % Time vector for response comparison
% t_compare = linspace(0, 5, 100);
% 
% % Define the input function u(t)
% u = @(t) sin(t);
% 
% % Step 1: Find Adjoint Matrix A*
% A_star = conj(A.');
% 
% % Step 2: Compute Lu_star(t) - Adjoint response for the adjoint system
% Lu_star = zeros(size(A_star, 1), length(t_compare));
% for i = 1:length(t_compare)
%     % Define the matrix exponential function
%     expAt_star = @(t) expm(A_star.' * t);
% 
%     % Define the system dynamics for the integral
%     dynamics_star = @(tau, z) expAt_star(t_compare(i) - tau) * (C.' * u(tau));
% 
%     % Initialize state variable for Lu_star
%     z0_star = zeros(size(A_star, 2), 1);
% 
%     % Integrate the system using a loop
%     for j = 1:length(t_compare)
%         [~, z_star] = ode45(@(tau, z) dynamics_star(tau, z), [eps t_compare(j)], z0_star); % Start from a small positive value
%         z0_star = z_star(end, :).';  % Update initial condition for the next step
%     end
% 
%     % Store the final state for Lu_star
%     Lu_star(:, i) = z0_star;
% end
% 
% % Compare Lu_star with Po
% figure;
% 
% subplot(2, 1, 1);
% plot(t_compare, Po, 'k', 'LineWidth', 2);
% title('Zero-Input Response (Po)');
% xlabel('Time');
% ylabel('State Variable');
% legend('Po');
% 
% subplot(2, 1, 2);
% plot(t_compare, Lu_star, 'k--', 'LineWidth', 2);
% title('Adjoint Response (Lu_star)');
% xlabel('Time');
% ylabel('State Variable');
% legend('Lu\_star');
% 
% % Display legend for comparison
% legend({'Po', 'Lu\_star'});