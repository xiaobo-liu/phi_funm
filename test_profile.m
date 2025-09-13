addpath('data','replication');
warning off

rng default
format compact 

nn = [20, 200, 500, 2500];
num_nn = length(nn);

pp = [1 4 7 10];
p = max(pp);

num_mat = 3;

Ttot = zeros(num_mat, num_nn);

Meval = zeros(num_mat, num_nn);
Mrecv = zeros(num_mat, num_nn);

Ppar = zeros(num_mat, num_nn);
Peval = zeros(num_mat, num_nn);
Precv = zeros(num_mat, num_nn);

main_loop = tic; % record the time consumption
for i = 1:num_mat

    for j = 1:num_nn

        n = nn(j);

        switch i
          case 1
            A = anymatrix('gallery/circul', n);		   % circulant matrix
          case 2
            A = anymatrix('gallery/triw', n, -2);  % upper triangular matrix
          case 3
            A = anymatrix('core/vand', n);         % Vandermonde matrix
        end

        [X, s, m, cost, runtime] = phi_funm_time(A, pp);
        
        Mrecv(i,j) = s*(p+1); % number of matmul in final recovering phase
        Meval(i,j) = round(cost - Mrecv(i,j) - 4/3); % number of matmul for the p+1 Pade approximents
        
        Ttot(i,j) = sum(runtime);
        
        Ppar(i,j) = runtime(1) / Ttot(i,j);
        Peval(i,j) = runtime(2) / Ttot(i,j);
        Precv(i,j) = runtime(3) / Ttot(i,j);
  
    end
end
filename = fullfile(pwd, 'data', 'test_profile.mat');
save(filename, 'nn', 'num_nn', 'pp', 'num_mat', 'Ttot', 'Meval', 'Mrecv', ...
    'Ppar', 'Peval', 'Precv');

fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);