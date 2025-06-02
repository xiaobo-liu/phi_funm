% A simple test 

addpath('replication')
format compact 

p = 10;
idx = 1:p;

% simple test 1
a = 2e10;
b = 4e8/6;
c = 200/3;
d = 3;
e = 1e-8;
A = [0 e 0; -(a+b) -d a; c 0 -c]; % Cleve Moler

fprintf('Running simple test 1...\n')
F = phi_funm(A,idx); 
X = phi_funm_ex(A,idx);

for i=1:length(idx)
    err_phi_funm = double(norm(X{i}-F{i},1)/norm(X{i},1));
    R = phipade(A,idx(i));
    err_phipade = double(norm(X{i}-R,1)/norm(X{i},1));
    fprintf('phi_funm_%d error: %.3e\n', idx(i), err_phi_funm);
    fprintf('phipade_%d error:  %.3e\n', idx(i), err_phipade);
end

% simple test 2
A = randn(5);
A = triu(A);
A(1,end) = 1e16;

fprintf('\n Running simple test 2...\n')
F = phi_funm(A,idx); 
X = phi_funm_ex(A,idx);

for i=1:length(idx)
    err_phi_funm = double(norm(X{i}-F{i},1)/norm(X{i},1));
    R = phipade(A,idx(i));
    err_phipade = double(norm(X{i}-R,1)/norm(X{i},1));
    fprintf('phi_funm_%d error: %.3e\n', idx(i), err_phi_funm);
    fprintf('phipade_%d error:  %.3e\n', idx(i), err_phipade);
end


function X = phi_funm_ex(A,idx)
% Compute reference solution using Symbolic Math Toolbox in 200 digits

num_digit = 200;
digits(num_digit);

p = max(idx);
n = size(A, 1);
I = speye(n);
E = zeros(n, n*p);
E(:, 1:n) = I;
J = tril(triu(ones(p), 1), 1);
J = kron(J, speye(n));
W = [A, E; zeros(n*p, n), J];

W = expm(vpa(W, num_digit));

% Output format
X = cell(1, numel(idx)); % Store all matrices in a cell array
for i = 1:numel(idx)
    j = idx(i);
    X{i} = W(1:n, j*n+1:(j+1)*n);
end

end