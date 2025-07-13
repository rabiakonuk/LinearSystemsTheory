% HW 8 
% Q1

clc
close all

% a
% Define the matrices A and B
A = [1 2 1; 0 4 3; 0 0 2];
B = [1; 1; 0];

% Number of states (rows of A)
n = size(A, 1);

% Initialize the controllability matrix C
C = B;

% Construct the controllability matrix
for i = 1:n-1
    C = [C, A^i * B];
end

% Check if the system is controllable by verifying the rank of C
isControllable = rank(C) == n;

% Display result
if isControllable
    disp('The system is controllable.');
else
    disp('The system is not controllable.');
end

% b
% Extract the first two columns of C, which span the controllable subspace
V = C(:, 1:2);

% Initialize the orthonormal basis vectors
orthonormalBasis = zeros(size(V));

% First orthonormal vector (normalized first vector of V)
orthonormalBasis(:, 1) = V(:, 1) / norm(V(:, 1));

% Second vector orthogonalized against the first
v2_orthogonal = V(:, 2) - proj(V(:, 2), orthonormalBasis(:, 1));

% Normalizing the second orthogonal vector
orthonormalBasis(:, 2) = v2_orthogonal / norm(v2_orthogonal);

% Display the orthonormal basis
disp('Orthonormal basis for the controllable subspace:');
disp(orthonormalBasis);

% c and d
% Compute e^(At) and integral_part (common for both parts c and d)
t = 1; % Time of transition
eAt = expm(A*t);
integral_part = integral(@(tau) expm(A*(t-tau)) * B, 0, t, 'ArrayValued', true);

% c
% Define the initial and final states for part c
x0_c = [1; 1; 1];
xf_c = [0; 0; 0];

% Solve for u(t) for part c
u_t_c = integral_part \ (xf_c - eAt*x0_c);

% Display the input for part c
disp('Input u(t) that drives the state from x(0) to x(1) for part c:');
disp(u_t_c);

% d
% Define the initial and final states for part d
x0_d = [0; 0; 0];
xf_d = [1; 1; 0];

% Solve for u(t) for part d
u_t_d = integral_part \ (xf_d - eAt*x0_d);

% Display the input for part d
disp('Input u(t) that drives the state from x(0) to x(1) for part d:');
disp(u_t_d);

% Function to calculate the projection of v onto u
function p = proj(v, u)
    p = (dot(v, u) / dot(u, u)) * u;
end