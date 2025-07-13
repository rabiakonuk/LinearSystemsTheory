% ME 564 HW9
% q5

% a
% Given system matrices
A = [0 1; -1 -2];
B = [1 0; 0 1];

% Desired eigenvalues
desired_eigenvalues = [0 0];

% Find state feedback matrix K using the place function
K = place(A, B, desired_eigenvalues);

% Display the state feedback matrix K
disp('State Feedback Matrix K:');
disp(K);

% Explanation:
% - desired_eigenvalues specifies the desired closed-loop evalues.
%   Here, we want both evalues to be at 0.
% - The place function calculates the state feedback matrix K s.t.
%   the closed-loop eigenvalues of A + BK match the desired_eigenvalues.
% - K is the state feedback matrix that satisfies the desired closed-loop
%   eigenvalues and is used in the state-feedback control law u(k) = Kx(k).

% b
% Extract the first row of B
B1 = B(1, :);

% Solve for the state feedback matrix K directly
K = acker(A, B1', desired_eigenvalues);

% Display the state feedback matrix K
disp('State Feedback Matrix K (for u1 only):');
disp(K);

% Explanation:
% - B1 is the first row of B, corresponding to the first input.
% - The acker function calculates the state feedback matrix K s.t.
%   the closed-loop eigenvalues of A + B1*K match the desired_eigenvalues.
% - K is the state feedback matrix for u1(k) = Kx(k) using only the first input.

% c
% Extract the second row of B
B2 = B(2, :);

% Solve for the state feedback matrix K directly for u2
K = acker(A, B2', desired_eigenvalues);

% Display the state feedback matrix K (for u2 only)
disp('State Feedback Matrix K (for u2 only):');
disp(K);

