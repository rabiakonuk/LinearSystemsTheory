% ME 564 HW7 Q5 c & d

% Define A(t)
A = @(t) [-4/t -2/t^2; 1 0];

% Define the state transition matrix function
% This is a numerical approximation, as MATLAB does not have a built-in
% symbolic matrix exponential for TVS
Phi = @(t, tau) expm(integral(@(s) A(s), tau, t, 'ArrayValued', true));

% Define the time vector (random)
tspan = [1 10]; 

% Define the initial condition
x0 = [1; 1];

% Solve the system using ode45
[t_out, x_out] = ode45(@(t, x) A(t)*x, tspan, x0);

% Display the solution at each time step
disp('Time and corresponding solution:')
disp([t_out, x_out])

% Compute and display the state transition matrix at a specific time
t = 10; 
tau = 1;
Phi_t_tau = Phi(t, tau);

disp(['State transition matrix Phi(' num2str(t) ', ' num2str(tau) '):'])
disp(Phi_t_tau)

% Compute and display the solution phi(t) using the state transition matrix
% This is the solution at time t starting from the initial condition x0 at time tau
phi_t = Phi_t_tau * x0;

disp(['Solution phi(' num2str(t) ') with initial condition x(1) = [1; 1]:'])
disp(phi_t)
