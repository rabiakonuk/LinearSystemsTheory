% ME564 HW9
% q4

% Given system matrices
A = [3 3 0 2; 0 87 0 60; 6 3 -3 2; 0 126 0 -87];
B = [0; 3; -1; 4];

% Controllability matrix
Co = ctrb(A, B);

% Check if the system is controllable
if rank(Co) == size(A, 1)
    % Perform Kalman decomposition
    Ac = A;
    Bc = B;
    Cc = eye(size(A));
    Tu = eye(size(A));
    Au = zeros(size(A));

    % Display the controllable and uncontrollable parts
    disp('Controllable Part:');
    disp(Ac);
    disp('Uncontrollable Part:');
    disp(Au);
else
    disp('System is not controllable');
end
