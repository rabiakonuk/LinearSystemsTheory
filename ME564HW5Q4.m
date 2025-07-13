% ME564_HW5
% Q4

% Define System Matrices
% The system is defined by the equation x(t+1) = Ax(t) + Bu(t) and y(t) = Cx(t)
% A is the state transition matrix, B is the input matrix, and C is the output matrix
A = [3, 0, -2; 0, 2, 5; 4, 3, -1];
B = [2, 0; 0, 0; 0, 1];
C = [1, 0, 1];

% Proposed Solution for Controllability
% The controllability of the system is checked using the controllability matrix.
% The controllability matrix is formed by [B, AB, A^2B, ..., A^(n-1)B]
% If the controllability matrix is of full rank, then the system is controllable.

% Calculate the controllability matrix
n = size(A, 1);  % Number of states
ControllabilityMatrix = [];
for i = 0:n-1
    ControllabilityMatrix = [ControllabilityMatrix, A^i * B];
end

% Check if the system is controllable
rank_C = rank(ControllabilityMatrix);
if rank_C == n
    disp('The system is controllable.');
else
    disp('The system is not controllable.');
end

% Proposed Solution for Observability
% The observability of the system is checked using the observability matrix.
% The observability matrix is formed by [C; CA; CA^2; ...; CA^(n-1)]
% If the observability matrix is of full rank, then the system is observable.

% Calculate the observability matrix
ObservabilityMatrix = [];
for i = 0:n-1
    ObservabilityMatrix = [ObservabilityMatrix; C * A^i];
end

% Check if the system is observable
rank_O = rank(ObservabilityMatrix);
if rank_O == n
    disp('The system is observable.');
else
    disp('The system is not observable.');
end
