% Compare the algebraic connectivity for the graph obtained with the optimal weights w* to the
% one obtained with w unif = (1/m)1 (i.e., a uniform allocation of weight to the edges).
%% data for maximizing algebraic connectivity of a graph
% provides A, F, g

rand('seed',0);
randn('seed',0);

n = 50; threshold = 0.2529;
rand('state',209);
xy = rand(n,2);

angle = 10*pi/180;
Rotate = [ cos(angle) sin(angle); -sin(angle) cos(angle) ];
xy = (Rotate*xy')';

Dist = zeros(n,n);
for i=1:(n-1);
  for j=i+1:n;
    Dist(i,j) = norm( xy(i,:) - xy(j,:) );
  end;
end;
Dist = Dist + Dist';
Ad = Dist < threshold;
Ad = Ad - eye(n);
m = sum(sum(Ad))/2;

% find the incidence matrix
A = zeros(n,m);
l = 0;
for i=1:(n-1);
  for j=i+1:n;
    if Ad(i,j)>0.5
      l = l + 1;
      A(i,l) =  1;
      A(j,l) = -1;
    end;
  end;
end;
A = sparse(A);

[n,m] = size(A);

F = ones(1,m);
g = 1;

% Matrix Q is a n x n-1 matrix and its columns are unit norm vectors, orthogonal to vector 1 and to each other
Q = null(ones(1,n));

cvx_begin
    % Weights
    variable w(m)
    % Laplacian of graph
    L = A*diag(w)*A';
    S = Q'*L*Q;
    % Objective function
    maximize (lambda_min(S))
    % Constraints
    subject to
    w >= 0;
    F*w <= g;
cvx_end
    
% Optimal weight vector v* has some zero entries (due to the finite precision of the solver => round small values of w* to exactly zero
w(abs(w) < 1e-4) = 0;

% Solution 1: eigenvalues with uniform allocation of weight to the edges
uniform_weights_solution = eig((1/m)*A*A');
% Solution 2: eigenvalues with uniform allocation of weight to the edges
optimal_weights_solution = eig(A*diag(w)*A');
% Algebraic connectivity = second smallest eigenvalue
eigen_value_index = 2;

fprintf('Algebraic connectivity with uniform constant weights: %f\n', uniform_weights_solution(eigen_value_index));
fprintf('Algebraic connectivity with optimal weights: %f\n', optimal_weights_solution(eigen_value_index));

% Plot graph with uniform weights
figure(1), clf % clear all graphics to create a new plot
gplot((1/m)*(A)*A', xy);
hold on;
plot(xy(:,1), xy(:,2), 'go','LineWidth',7, 'MarkerSize',3);
axis([0.05 1.1 -0.05 1.05]);
% title('Graph topology with uniform weights')
hold off;

% Plot graph with optimal weights
figure(2), clf % clear all graphics to create a new plot
gplot(A*diag(w)*A',xy);
hold on;
plot(xy(:,1), xy(:,2), 'go','LineWidth',7, 'MarkerSize',3);
axis([0.05 1.1 -0.05 1.05]);
% title('Graph topology with optimal weights')
hold off;

% Plot graph with edges thickness representing weights - optimal
figure(3), clf
plotgraph(A,xy,w);

% Plot graph with edges thickness representing weights - uniform
unif_weights = (1/m) * ones(m);
figure(4), clf
plotgraph(A,xy,unif_weights);
