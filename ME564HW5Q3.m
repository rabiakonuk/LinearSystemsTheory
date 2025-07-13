% ME564_HW5
% Q3

% Part a
disp('Part a:');

% i. Find the eigenvalues
B = [8, -8, -2; 4, -3, -2; 3, -4, 1];
eigenvalues = eig(B);
disp('Eigenvalues:');
disp(eigenvalues);

% ii. Find eigenvectors and/or generalized eigenvectors
[V, J] = jordan(B);
disp('Eigenvectors/Generalized Eigenvectors (columns of P):');
disp(V);

% iii. Compute the Jordan form
J_computed = inv(V) * B * V;
disp('Jordan Form:');
disp(J_computed);

% Double-check with MATLAB's jordan function
[J_check, P_check] = jordan(B);
disp('Jordan Form (MATLAB check):');
disp(J_check);

% Part b
disp('Part b:');

% i. Find the eigenvalues
B = [2, 1, 1; 0, 3, 1; 0, -1, 1];
eigenvalues = eig(B);
disp('Eigenvalues:');
disp(eigenvalues);

% ii. Find eigenvectors and/or generalized eigenvectors
[V, J] = jordan(B);
disp('Eigenvectors/Generalized Eigenvectors (columns of P):');
disp(V);

% iii. Compute the Jordan form
J_computed = inv(V) * B * V;
disp('Jordan Form:');
disp(J_computed);

% Double-check with MATLAB's jordan function
[J_check, P_check] = jordan(B);
disp('Jordan Form (MATLAB check):');
disp(J_check);
