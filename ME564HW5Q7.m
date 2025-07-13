% ME564_HW5
% Q7 

% a
% Define the matrix A
A = [2, 1; 1, -8];

% Initialize the approximation for e^A as the identity matrix
approx_eA = eye(size(A));

% Compute the approximation using the first five terms of the series definition
terms = 5;
for k = 1:terms
    approx_eA = approx_eA + (A^k) / factorial(k);
end

disp('Approximation for e^A using the first five terms:');
disp(approx_eA);

% b
% Initialize the approximation for e^-A as the identity matrix
approx_e_neg_A = eye(size(A));

% Compute the approximation using the first five terms of the series definition
for k = 1:terms
    approx_e_neg_A = approx_e_neg_A + ((-1)^k) * (A^k) / factorial(k);
end

% Compute the approximation for e^A as the inverse of approx_e_neg_A
approx_eA = inv(approx_e_neg_A);

disp('Approximation for e^A using the first five terms and inverse method:');
disp(approx_eA);

% c
% Compute the eigenvalues and eigenvectors
[V, D] = eig(A);

disp('Eigenvalues of A:');
disp(diag(D));
disp('Eigenvectors of A (columns):');
disp(V);

% d
% Compute e^D
eD = exp(diag(D));

% Compute e^A using the formula e^A = P e^D P^{-1}
eA_exact = V * diag(eD) / V;

% Compute e^A using MATLAB's expm function
eA_matlab = expm(A);

% Display the results
disp('Exact value of e^A:');
disp(eA_exact);
disp('Value of e^A using MATLAB''s expm function:');
disp(eA_matlab);

% e
% Construct matrix P using the eigenvectors
P = V;

% Construct diagonal matrix Λ using the eigenvalues
Lambda = D;

% Compute e^Λ
eLambda = exp(diag(Lambda));

% Compute e^A using the formula e^A = P * (e^Λ) * P^(-1)
eA_diagonalization = P * diag(eLambda) / P;

% Display the results
disp('Value of e^A using diagonalization:');
disp(eA_diagonalization);
disp('Value of e^A using MATLAB''s expm function for comparison:');
disp(expm(A));