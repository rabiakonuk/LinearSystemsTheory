% HW 8
% Q3

% Define symbolic variables for the system
syms a11 a12 a13 a21 a22 a23 a31 a32 a33
syms b1 b2 b3
syms x01 x02 x03 xf1 xf2 xf3

% Define the matrices A, B, x0 and xf using symbolic variables
A = [a11, a12, a13; a21, a22, a23; a31, a32, a33];
B = [b1; b2; b3];
x0 = [x01; x02; x03];
xf = [xf1; xf2; xf3];

% Number of states (rows of A)
n = size(A, 1);

% Construct the controllability matrix C
C = B;
for i = 1:n-1
    C = [C, A^i * B];
end

% Assuming the system is controllable and C is full rank
% the state transition from x0 to xf
x_transition = A^(n-1) * x0;

% the input sequence
u = C \ (xf - x_transition); % Using the matrix division operator

% Convert symbolic expressions to strings and display the first 10 terms
disp('Input sequence u that drives the state from x0 to xf (first terms):');
for i = 1:min(10, length(u))
    disp(['u(', num2str(i-1), ') = ', char(u(i))]);
end
if length(u) > 2
    disp('...');
end
