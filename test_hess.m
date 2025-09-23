% Test on Hessenberg matrices.

% Create directories to store the results, if not exist
if ~exist('figs', 'dir'), mkdir('figs'); end
if ~exist('data', 'dir'), mkdir('data'); end

addpath('data','replication','figs')

rng default

mats_id = {'bcspwr10', 'gr_30_30', 'helm2d03', 'orani678', 'poisson99'};
num_mat = 5;

mm = [30, 80];
num_mm = length(mm);

pp = [1 4];
num_pp = length(pp);


cost_phifunc = zeros(num_mat, num_mm, num_pp);
cost_phipade_dft = zeros(num_mat, num_mm, num_pp);
cost_phipade_opt = zeros(num_mat, num_mm, num_pp);

error_phifunc = zeros(num_mat, num_mm, num_pp);
error_phipade_dft = zeros(num_mat, num_mm, num_pp);
error_phipade_opt = zeros(num_mat, num_mm, num_pp);

main_loop = tic; % record the time consumption
for i = 1:num_mat
    switch i
        case 1
            mat = load('bcspwr10');
            A = mat.Problem.A;
            
        case 2
            mat = load('gr_30_30');
            A = mat.Problem.A;

        case 3
            mat = load('helm2d03');
            A = mat.Problem.A;

        case 4
            mat = load('orani678');
            A = mat.Problem.A;

        case 5
            A = - 2500 * anymatrix('gallery/poisson', 99);

    end
    
    n = size(A, 1);
    b = ones(n, 1);

    for j = 1:num_mm
        m = mm(j);
        [Q, H] = arnoldi(A, b, m);
        % I = eye(m); 
        % res = A*Q(:,1:m) - Q(:,1:m)*H(1:m,1:m) - H(m+1,m)*Q(:,m+1)*I(:,m)';
        % norm(res)

        H = H(1:m, 1:m);
        
        X_phipade_dft = zeros(num_pp, m, m);
        X_phipade_opt = zeros(num_pp, m, m);
        X_phifunc = zeros(num_pp, m, m);
        
        X = phi_func_ex(H, pp); % reference solution

        for k = 1:num_pp
            fprintf('Running the test...Matrix id: %s\n m = %d, p = %d\n', ...
                mats_id{i}, mm(j), pp(k));

            p = pp(k);

            [X_phifunc_cell, ~, ~, cost_phifunc(i,j,k)] = phi_funm(H, p);
            X_phifunc(k,:,:) = cell2mat(X_phifunc_cell);
         
            dft_deg_phipade = 7; % phipade default degree is m = 7

            if p==1
                
                % phipade degree and scaling optimized for cost
                [deg_phipade, cost_phipade_opt(i,j,k)] = select_deg_phipade(H, p);
                G1 = phipade(H, p, deg_phipade); 
                X_phipade_opt(k,:,:) = G1;
  
             
                F1 = phipade(H, p, dft_deg_phipade); 
                cost_phipade_dft(i,j,k) = phipade_default_cost(H, p);

                X_phipade_dft(k,:,:) = F1;

            else % p = 4

                % phipade degree and scaling optimized for cost
                [deg_phipade, cost_phipade_opt(i,j,k)] = select_deg_phipade(H, p);
                [G1, G2, G3, G4] = phipade(H, p, deg_phipade); 
                X_phipade_opt(k,:,:) = G4;
 
                [F1, F2, F3, F4] = phipade(H, p, dft_deg_phipade); 
                cost_phipade_dft(i,j,k) = phipade_default_cost(H, p);
                
                X_phipade_dft(k,:,:) = F4;
   
            end
            
            normXk = norm(X{k},1);
            error_phipade_opt(i,j,k) = double( norm( X{k}-squeeze(X_phipade_opt(k, :, :) ), 1) / normXk );
            error_phipade_dft(i,j,k) = double( norm( X{k}-squeeze(X_phipade_dft(k, :, :) ), 1) / normXk );
            error_phifunc(i,j,k) = double( norm( X{k}-squeeze(X_phifunc(k,:,:) ), 1) / normXk );
  
        end
    end
end

filename = fullfile(pwd, 'data', 'test_hess.mat');
save(filename, 'mats_id', 'num_mat', 'mm', 'num_mm', 'pp', 'num_pp', ...
    'cost_phifunc', 'cost_phipade_dft', 'cost_phipade_opt', ...
    'error_phifunc', 'error_phipade_dft', 'error_phipade_opt');

fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);