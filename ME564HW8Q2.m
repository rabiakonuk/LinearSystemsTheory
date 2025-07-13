% HW 8 
% Q2

clc
close all

% a
% Define the matrices A and B
A = [-1 0 -1; 0 -3 1; 0 0 -2];
B = [1; 0; 2];

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
% Extract the columns of C that span the controllable subspace
V = C(:, 1:rank(C));

% Initialize the orthonormal basis vectors
orthonormalBasis = zeros(size(V));

% First orthonormal vector (normalized first vector of V)
orthonormalBasis(:, 1) = V(:, 1) / norm(V(:, 1));

% Second vector orthogonalized against the first and normalized
v2_orthogonal = V(:, 2) - proj(V(:, 2), orthonormalBasis(:, 1));
orthonormalBasis(:, 2) = v2_orthogonal / norm(v2_orthogonal);

% Third vector orthogonalized against the first two and normalized
v3_orthogonal = V(:, 3) - proj(V(:, 3), orthonormalBasis(:, 1)) - proj(V(:, 3), orthonormalBasis(:, 2));
orthonormalBasis(:, 3) = v3_orthogonal / norm(v3_orthogonal);

% Display the orthonormal basis
disp('Orthonormal basis for the controllable subspace:');
disp(orthonormalBasis);

% c
% Define the time interval
t0 = 0;
t1 = 1;

% Calculate the controllability Gramian
Wc = integral(@(tau) expm(A*(t1-tau))*B*B'*expm(A'*(t1-tau)), t0, t1, 'ArrayValued', true);

% Display the controllability Gramian
disp('Controllability Gramian from t = 0 to t = 1:');
disp(Wc);

% d
% Define the initial and final states
x0 = [1; 1; 1];
xf = [0; 0; 0];

% Time of transition
t = 1;

% Compute e^(At)
eAt = expm(A*t);

% Solving for the integral part (as in previous part)
integral_part = integral(@(tau) expm(A*(t-tau)) * B, 0, t, 'ArrayValued', true);

% Solve for u(t)
u_t = integral_part \ (xf - eAt*x0);

% Display the input
disp('Input u(t) that drives the state from x(0) to x(1):');
disp(u_t);

% e
% Define Q as B*B'
Q = B * B';

% Solve the Lyapunov equation A'P + PA + Q = 0 for P
P = lyap(A, Q);

% Display the P matrix
disp('Matrix P that solves the Lyapunov equation:');
disp(P);

% f
% Compute the infinite-time controllability Gramian (which should be equal to P)
% AW + WA' + BBt = 0
Wc_infinite = lyap(A, B*B');

% Display the infinite-time controllability Gramian
disp('Controllability Gramian from t = 0 to t = infinity:');
disp(Wc_infinite);


% Function to calculate the projection of v onto u
function p = proj(v, u)
    p = (dot(v, u) / dot(u, u)) * u;
end

