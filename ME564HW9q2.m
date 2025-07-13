% ME564 HW9
% q2

% Transfer function G(s)
numerator = [1 3];
denominator = conv([1 2], [1 1]);
G = tf(numerator, denominator);

% Function to display matrices
displayMatrices = @(A, B, C, D, label) fprintf('%s:\nMatrix A:\n%s\n\nMatrix B:\n%s\n\nMatrix C:\n%s\n\nMatrix D:\n%s\n\n', label, mat2str(A), mat2str(B), mat2str(C), mat2str(D));

% State-space representation
[A, B, C, D] = ssdata(ss(G));

% a - controllable canonical realization
displayMatrices(A, B, C, D, 'Controllable Canonical Realization');

% b - observable canonical realization
n = size(A, 1); % Number of states
Ao = [A; eye(n-1) zeros(n-1, 1)]; % Extend A matrix
Bo = [zeros(1, n-1) 1]'; % Extend B matrix
Co = C; % C matrix remains the same
displayMatrices(Ao, Bo, Co, D, 'Observable Canonical Realization');

% c - modal realization
Am = diag(eig(A)); % Diagonal matrix with eigenvalues of A
Bm = B; % B matrix remains the same
Cm = C; % C matrix remains the same
displayMatrices(Am, Bm, Cm, D, 'Modal Realization');

% d - non-minimal realization
Anm = [A, zeros(size(A, 1), 1); zeros(1, size(A, 2)), -1]; % Extend A to 4x4
Bnm = [B; 0]; % Extend B to 4x1
Cnm = [C, zeros(1, size(C, 2))]; % Extend C to 1x4
displayMatrices(Anm, Bnm, Cnm, D, 'Non-minimal Realization');
