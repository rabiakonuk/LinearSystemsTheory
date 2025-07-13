% ME564_HW5
% Q1

% a
% Part i: Find Eigenvalues and Eigenvectors
A = [2, 3; 3, 1];
[V, D] = eig(A);

% Part ii: Plot the Unit Ball and Its Transformation
theta = linspace(0, 2*pi, 100);
x = cos(theta);
y = sin(theta);
unitBall = [x; y];

transformedBall = A * unitBall;

figure;
subplot(1,2,1);
plot(x, y);
title('Unit Ball for Matrix A');
axis equal;

subplot(1,2,2);
plot(transformedBall(1,:), transformedBall(2,:));
title('Transformed Unit Ball for Matrix A');
axis equal;

% Part iii: Plot the Eigenvectors
maxLength = 5;
length = linspace(0, maxLength, 100);

figure;
hold on;
for i = 1:size(V, 2)
    v = V(:, i);
    plot(length * v(1), length * v(2));
end
title('Eigenvectors (a)');
axis equal;
hold off;

% Part iv: Find the Value of maxLength
eigenvalues = diag(D);
maxLengthValues = 1 ./ eigenvalues;

disp('Values of maxLength for each eigenvector (a):');
disp(maxLengthValues);

% b

% answer to d: The reason eigenvectors were not plotted for part (b) is 
% likely because the matrix A in that part is a rotation matrix & for a 
% 2D rotation matrix, the eigenvalues & eigenvectors are complex. 
% So the reason might be: Plotting complex eigenvectors in the same 
% 2D space as the unit ball and its transformation would not be meaningful, 
% as the eigenvectors would not lie in the same real plane.

% Part i: Find Eigenvalues and Eigenvectors for Matrix B
A = [cos(pi/5), -sin(pi/5); sin(pi/5), cos(pi/5)];
[V, D] = eig(A);

% Display Eigenvalues and Eigenvectors
disp('Eigenvalues (b):');
disp(diag(D));
disp('Eigenvectors (b):');
disp(V);

% Part ii: Plot the Unit Ball and Its Transformation for Matrix B
theta = linspace(0, 2*pi, 100);
x = cos(theta);
y = sin(theta);
unitBall = [x; y];

transformedBall = A * unitBall;

figure;
subplot(1,2,1);
plot(x, y);
title('Unit Ball for Matrix B');
axis equal;

subplot(1,2,2);
plot(transformedBall(1,:), transformedBall(2,:));
title('Transformed Unit Ball for Matrix B');
axis equal;

% c
% Part i: Find Eigenvalues and Eigenvectors for Matrix C
A = [7/8, -1/4; -1/8, 1];
[V, D] = eig(A);

% Display Eigenvalues and Eigenvectors
disp('Eigenvalues (c):');
disp(diag(D));
disp('Eigenvectors (c):');
disp(V);

% Part ii: Plot the Unit Ball and Its Transformation for Matrix C
theta = linspace(0, 2*pi, 100);
x = cos(theta);
y = sin(theta);
unitBall = [x; y];

transformedBall = A * unitBall;

figure;
subplot(1,2,1);
plot(x, y);
title('Unit Ball for Matrix C');
axis equal;

subplot(1,2,2);
plot(transformedBall(1,:), transformedBall(2,:));
title('Transformed Unit Ball for Matrix C');
axis equal;

% Part iii: Plot the Eigenvectors for Matrix C
maxLength = 5;
length = linspace(0, maxLength, 100);

figure;
hold on;
for i = 1:size(V, 2)
    v = V(:, i);
    plot(length * v(1), length * v(2));
end
title('Eigenvectors for Matrix C');
axis equal;
hold off;

% Part iv: Find the Value of maxLength for Matrix C
eigenvalues = diag(D);
maxLengthValues = 1 ./ eigenvalues;

disp('Values of maxLength for each eigenvector for Matrix C:');
disp(maxLengthValues);
