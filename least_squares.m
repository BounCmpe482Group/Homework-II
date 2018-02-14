clc;clear;

% load data
load hw2_data.mat
[n,~] = size(X);

% Formulate the given problem as a least squares problem. (Given A, y
% find x such that (Ax-y) is minimized.)
A = [ones(n-3,1) X(3:n-1,2) X(2:n-2,2) X(1:n-3,2)];
y = X(4:end,2);

% Find QR decomposition, store it in A (in-place)
A = qr_householder(A);

[~, n_A] = size(A);

% Calculate Q'*b. Apply the matrices Q_1, Q_2, ...
% Q_n to y to obtain x (again, in-place)
for k = 1:n_A
    v_k = A(k+1:end,k);
    y(k:end,1) = y(k:end,1) - 2*v_k*(v_k'*y(k:end,1));
end

% We only need first n_A elements of y
y = y(1:n_A, 1);

% Now that we have Q'*b and we also know Rx = Q'*b, multiply Q'*b with R
% inverse to get x (or in the context of the problem, beta)

%Get upper triangular matrix R_hat
R = triu(A);
R = R(1:n_A, 1:n_A);

% Multiply with inverse
beta = inv(R)*y;

% Print result
disp(beta);

% Un-comment to write the results to file
%fileID = fopen('estimates.txt','wt');
%fprintf(fileID,'%f\n',beta);