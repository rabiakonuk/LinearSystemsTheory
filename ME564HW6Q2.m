% ME564 HW6 Q1
clc;
clear all;

% a
% Given matrix A and vector b
A = [3, 0, 2, -2; 0, 2, -1, -1; 6, 0, 4, -4];
b = [1; 1; 1];

% Compute A* * A
A_star_A = A' * A;

% Perform SVD on A* * A
[U, S, V] = svd(A_star_A);

% Extract non-zero singular values and corresponding vectors to form V1
non_zero_indices = find(diag(S) > 1e-10);
V1 = V(:, non_zero_indices);

% Display V1
disp('Orthonormal basis V1 for R(A*):');
disp(V1);

% b
% Perform SVD on A
[U, S, V] = svd(A);

% Extract zero singular values and corresponding vectors to form V2
zero_indices = find(diag(S) <= 1e-10);
V2 = V(:, zero_indices);

% Display V2
disp('Orthonormal basis V2 for N(A):');
disp(V2);

% c
% Extract singular values
singular_values = diag(S);

% Display singular values
disp('Singular values of A:');
disp(singular_values);

% d
% Extract non-zero singular values and corresponding vectors to form V1
non_zero_indices = find(diag(S) > 1e-10);
V1 = V(:, non_zero_indices);
S_non_zero = diag(S(non_zero_indices, non_zero_indices));

% Compute S^-1
S_inv = diag(1 ./ S_non_zero);

% Compute U1 = A * V1 * S^-1
U1 = A * V1 * S_inv;

% Verify orthonormality
orthonormal_check = U1' * U1;

% Display U1 and orthonormality check
disp('Computed U1:');
disp(U1);
disp('Orthonormality check (U1'' * U1):');
disp(orthonormal_check);

% e
% Initialize U2
U2 = [];

% Use Gram-Schmidt to find vectors orthogonal to U1
for i = 1:3
    u = eye(3, i);
    for j = 1:size(U1, 2)
        u = u - (U1(:, j)' * u) .* U1(:, j);
    end
    u = u / norm(u);
    U2 = [U2, u];
end

% Display U2
disp('Orthonormal basis U2:');
disp(U2);

% f
% Compute the pseudo-inverse A_dagger = V1 * S_inv * U1'
A_dagger = V1 * S_inv * U1';

% Compute the least squares solution x = A_dagger * b
x_ls = A_dagger * b;

% Display the least squares solution
disp('Least squares solution x:');
disp(x_ls);