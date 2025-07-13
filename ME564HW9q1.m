% ME564 - HW9
% Q1

% Define system matrices
A = [2 1 1; 5 3 6; -5 -1 -4];
B = [1; 0; 0];
C = [1 1 2];

% (a) Identify the three modes of the system

% Calculate eigenvalues and eigenvectors
[eigVectors, eigValues] = eig(A);

% Extract the eigenvalues into a vector
eigValues = diag(eigValues);

% Display the eigenvalues and corresponding eigenvectors
disp('Eigenvalues:');
disp(eigValues);
disp('Eigenvectors:');
disp(eigVectors);

% Check for defective eigenvalues and display modes
for i = 1:length(eigValues)
    algebraicMultiplicity = sum(eigValues(i) == eigValues);
    geometricMultiplicity = rank(eigVectors(:, eigValues == eigValues(i)));
    
    if algebraicMultiplicity == geometricMultiplicity
        % Non-defective eigenvalue
        mode = eigVectors(:, i) * exp(eigValues(i) * 't');
        disp(['Mode for eigenvalue ', num2str(eigValues(i)), ':']);
        disp(mode);
    else
        % Defective eigenvalue, attempt to find generalized eigenvector
        genEigVector = null((A - eigValues(i) * eye(size(A)))^2, 'r');
        if ~isempty(genEigVector)
            % Generalized eigenvector found
            mode = genEigVector * exp(eigValues(i) * 't');
            disp(['Defective mode for eigenvalue ', num2str(eigValues(i)), ':']);
            disp(mode);
        else
            % Unable to find a generalized eigenvector
            disp(['Eigenvalue ', num2str(eigValues(i)), ' is defective, but unable to find a generalized eigenvector.']);
        end
    end
end

% (b) Identify whether each mode is controllable or uncontrollable

% Calculate the controllability matrix
n = size(A, 1);
controllabilityMatrix = B;
for i = 1:n-1
    controllabilityMatrix = [controllabilityMatrix, A^i * B];
end

% Check controllability for each mode
for i = 1:n
    eigVec = eigVectors(:, i);
    if rank([controllabilityMatrix, eigVec]) == rank(controllabilityMatrix)
        disp(['Mode corresponding to eigenvector ', num2str(i), ' is controllable.']);
    else
        disp(['Mode corresponding to eigenvector ', num2str(i), ' is uncontrollable.']);
    end
end

% (c) Identify whether each mode is observable or unobservable

% Calculate the observability matrix
observabilityMatrix = C;
for i = 1:n-1
    observabilityMatrix = [observabilityMatrix; C * A^i];
end

% Check observability for each mode
for i = 1:n
    eigVec = eigVectors(:, i);
    if rank([observabilityMatrix; eigVec']) == rank(observabilityMatrix)
        disp(['Mode corresponding to eigenvector ', num2str(i), ' is observable.']);
    else
        disp(['Mode corresponding to eigenvector ', num2str(i), ' is unobservable.']);
    end
end

% (d) Transfer Function Y(s)/U(s) Confirmation

% Define the symbolic variable s
syms s

% Compute the transfer function
transferFunction = C * inv(s * eye(size(A)) - A) * B;

% Simplify and display the transfer function
transferFunction = simplify(transferFunction);
disp('Transfer Function Y(s)/U(s):');
disp(transferFunction);

% (e) Transform the system into modal form

% Find the Jordan normal form J of matrix A and the change of basis matrix V
[V, J] = jordan(A);

% Compute the transformed B and C matrices
B_bar = inv(V) * B;
C_bar = C * V;

% Display the results
disp('Jordan normal form J:');
disp(J);
disp('Change of basis matrix V:');
disp(V);
disp('Transformed B matrix (B_bar):');
disp(B_bar);
disp('Transformed C matrix (C_bar):');
disp(C_bar);
