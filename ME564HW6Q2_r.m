% ME564 HW6 Q2

clc;
clear all;

% Given matrix A
A = [0.0974, -0.1178, 0.7876, -0.1168, 0.0178;
     0.1291, -0.1174, 1.2850, 0.0302, 0.0971;
     0.0528, 0.1119, 0.1325, 0.7668, 0.0637;
     0.0424, 0.2647, 0.2806, 1.7644, 0.1195];

% a
% Perform SVD on A
[U, S, V] = svd(A);

% Display U, S, and V
disp('Matrix U:');
disp(U);
disp('Matrix S (Sigma):');
disp(S);
disp('Matrix V:');
disp(V);

% Extract and display singular values
singular_values = diag(S);
disp('Singular values:');
disp(singular_values);

% b
% Extract and sort singular values
[sorted_singular_values, indices] = sort(singular_values, 'descend');

% Select the two largest singular values to form Sigma_r
Sigma_r = diag(sorted_singular_values(1:2));

% Extract corresponding columns from U and V to form U_r and V_r
U_r = U(:, indices(1:2));
V_r = V(:, indices(1:2));

% Display U_r, Sigma_r, and V_r
disp('Reduced Matrix U_r:');
disp(U_r);
disp('Reduced Matrix Sigma_r:');
disp(Sigma_r);
disp('Reduced Matrix V_r:');
disp(V_r);

% c
% Compute A_r = U_r * Sigma_r * V_r'
A_r = U_r * Sigma_r * V_r';

% Display the reduced-order matrix A_r
disp('Reduced-order matrix A_r:');
disp(A_r);

% d
% Compute the difference between A and A_r
A_diff = A - A_r;

% Compute the Frobenius norm of the difference
Frobenius_norm = norm(A_diff, 'fro');

% Display the Frobenius norm
disp('Frobenius norm of (A - A_r):');
disp(Frobenius_norm);