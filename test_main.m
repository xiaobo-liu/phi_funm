% Test against phipade and other algorithms.

addpath('replication')
anymatrix scan
rng default
format compact 

pp = [1 4 7 10]; % here for p > 0  
pmax = max(pp);

num_pp = length(pp);

u = eps('double')/2;
ids_min = 1;
[~, ids_max] = testmats();
num_mats = ids_max - ids_min + 1;

err_phifunc = zeros(num_pp, num_mats);
err_phipade_opt = zeros(num_pp, num_mats); % phipade with optimized degree
err_phipade_dft = zeros(num_pp, num_mats); % phipade with default degree
err_expmblktri = zeros(num_pp, num_mats);

% for p = 0
err_expm = zeros(1, num_mats);
err_phifunc0 = zeros(1, num_mats);

cost_phifunc = zeros(num_mats,1);
cost_phipade_opt = zeros(num_mats,1);
cost_phipade_dft = zeros(num_mats,1);

time_phifunc = zeros(num_mats,1);
time_phipade_opt = zeros(num_mats,1);
time_phipade_dft = zeros(num_mats,1);

condest1u_expm = zeros(1, num_mats);
condest1u = zeros(num_pp, num_mats);
        
main_loop = tic; % record the time consumption
for k = ids_min:ids_max
    n = 20; % default matrix size
    fprintf('Running the test...Matrix id: %d\n', k);
    A = testmats(k,n);

    n = size(A,1);
    I = eye(n);
    X = phi_func_ex(A, pp, 0); % reference solution, X{end} for p = 0
    tic;
    [X_phifunc, ~, ~, cost_phifunc(k)] = phi_funm(A, pp, 0); % X_phifunc{end} for p=0
    time_phifunc(k) = toc;

    % phipade degree and scaling optimized for cost 
    [deg_phipade, cost_phipade_opt(k)] = select_deg_phipade(A, pmax);
    tic;
    [G1, G2, G3, G4, G5, G6, G7, G8, G9, G10] = phipade(A, pmax, deg_phipade); 
    time_phipade_opt(k) = toc;

    X_phipade_opt = zeros(num_pp, n, n);
    X_phipade_opt(1,:,:) = G1;
    X_phipade_opt(2,:,:) = G4;
    X_phipade_opt(3,:,:) = G7;
    X_phipade_opt(4,:,:) = G10;

    % phipade default degree is m = 7
    dft_deg_phipade = 7;
    tic;
    [F1, F2, F3, F4, F5, F6, F7, F8, F9, F10] = phipade(A, pmax, dft_deg_phipade); 
    time_phipade_dft(k) = toc;

    cost_phipade_dft(k) = phipade_default_cost(A, pmax);
    X_phipade_dft = zeros(num_pp, n, n);
    X_phipade_dft(1,:,:) = F1;
    X_phipade_dft(2,:,:) = F4;
    X_phipade_dft(3,:,:) = F7;
    X_phipade_dft(4,:,:) = F10;

    X_expm = expm(A); 

    % for expm_blktri, W = [A E; zeros(pmax*n, n); J];
    E = [I zeros(n, (pmax-1)*n)];
    Jp = gallery('jordbloc', pmax, 0);
    J = kron(Jp, I);
    [~, ~, X_expmblktri_full] = expm_block_tri(A, J, E); 
        
    % for p = 0
    expm_ex = X{end};
    err_phifunc0(k) = double( norm(expm_ex-X_phifunc{end}, 1) / norm(expm_ex, 1) );
    err_expm(k) = double( norm(expm_ex-X_expm, 1) / norm(expm_ex, 1) );
    funm = @(A,E)funm_fd(@expm,A,E);
    condest1u_expm(k) = funm_condest1(A, @expm,funm) * u;
        
    % for p > 0
    for j=1:num_pp
        X_expmblktri = X_expmblktri_full(:,(pp(j)-1)*n+1:pp(j)*n);         
        normXj = norm(X{j},1);
        err_phipade_opt(j, k) = double( norm(X{j}-squeeze(X_phipade_opt(j, :, :)), 1) / normXj );
        err_phipade_dft(j, k) = double( norm(X{j}-squeeze(X_phipade_dft(j, :, :)), 1) / normXj );
        err_phifunc(j, k) = double( norm(X{j}-X_phifunc{j}, 1) / normXj );
        err_expmblktri(j, k) = double( norm(X{j}-X_expmblktri, 1) / normXj );
        phi_funcj = @(A) cell2mat(phi_funm(A, pp(j)));

        funm = @(A,E)funm_fd(phi_funcj,A,E); % Awad
        % Computing FD via funm_fd is more accurate than finite diff.
        % used by funm_condest1 if nargin < 3.
        condest1u(j, k) = funm_condest1(A, phi_funcj, funm) * u; % Awad
    end     
end
fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);
    
%% plot macros

clf;
lg_lindwidth = 1.8;
lg_markersize = 5;
lg_fontsize = 14;
title_fontsize = 14;

axlabel_lindwidth = 1.0;
axlabel_fontsize = 10;

color_cond    = [0 0 0];
color_expm = [0.635 0.078 0.184];
color_phi_funm = [0.23 0.48 0.34]; 
color_expmblktri = [0.635 0.078 0.184];
color_phipade_dft  = [0.8500, 0.3250, 0.098];
color_phipade_opt  = [0, 0.4470, 0.7410];

%% for p = 0
[~, perm] = sort(condest1u_expm, 'descend');

% relative errors
figure(1)
ga = gobjects(3, 1);

ga(1) = semilogy(1:num_mats, condest1u_expm(perm), '-', ...
    'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize, 'DisplayName', '$\kappa_{\exp}(A)u$');
hold on 
ga(2) = semilogy(1:num_mats, err_expm(perm), '^', ...
       'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize, 'DisplayName', '\texttt{expm}');
hold on
ga(3) = semilogy(1:num_mats, err_phifunc0(perm), 'v', ...
    'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize, 'DisplayName', '\texttt{phi\_funm}');

mycolors = [color_cond; color_expm; color_phi_funm];
ax = gca; 
ax.ColorOrder = mycolors;
legend(ga([1,2,3]), 'NumColumns', 1, 'FontSize', lg_fontsize, 'interpreter', 'latex');

grid on;

set(legend, 'Location', 'NorthEast');

set(gca,'linewidth', axlabel_lindwidth)
set(gca,'fontsize', axlabel_fontsize)
xlim([1,num_mats]);
xticks([1 20:20:100 num_mats]);
    
ylim([1e-18 1])
yticks(10.^(-18:3:0))

% performance profile
figure(2)
ratio_max = 12;   
markers = {'-v', '-^'};
err_matix = [err_phifunc0' err_expm'];
methods_names = {'\texttt{phi\_funm}', '\texttt{expm}'};

myperfprof(err_matix, ratio_max, flip(mycolors(2:end,:),1), markers, lg_markersize, lg_lindwidth);
legend(methods_names, 'interpreter', 'latex', 'FontSize', lg_fontsize);
grid on;
ylim([0.2 1])

set(legend, 'Location', 'southeast');
set(gca,'linewidth', axlabel_lindwidth)
set(gca,'fontsize', axlabel_fontsize)


%% p > 0
for j=1:num_pp
    err_phifunc_j = err_phifunc(j,:);
    err_phipade_opt_j = err_phipade_opt(j,:);
    err_phipade_dft_j = err_phipade_dft(j,:);
    err_expmblktri_j = err_expmblktri(j,:);
    condest1u_j = condest1u(j, :);

    [~, perm] = sort(condest1u_j, 'descend');
   
    % relative errors
    figure(2*j+1)  

    ga = gobjects(5, 1);

    ga(1) = semilogy(1:num_mats, condest1u_j(perm), '-', ...
        'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize, 'DisplayName', '$\kappa_{\varphi_j}(A)u$');
    hold on 
    ga(2) = semilogy(1:num_mats, err_expmblktri_j(perm), '^', ...
        'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize, 'DisplayName', '\texttt{expm\_blktri}');
    hold on 
    ga(3) = semilogy(1:num_mats, err_phipade_dft_j(perm), 's', ...
        'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize, 'DisplayName', '\texttt{phipade\_dft}');
    hold on
    ga(4) = semilogy(1:num_mats, err_phipade_opt_j(perm), 'o', ...
        'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize, 'DisplayName', '\texttt{phipade\_opt}');
    hold on
    ga(5) = semilogy(1:num_mats, err_phifunc_j(perm), 'v', ...
        'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize, 'DisplayName', '\texttt{phi\_funm}');
    mycolors = [color_cond; color_expmblktri; color_phipade_dft; color_phipade_opt; color_phi_funm];

    ax = gca; 
    ax.ColorOrder = mycolors;
    legend(ga([1,2,3,4,5]), 'NumColumns', 1, 'FontSize', lg_fontsize, 'interpreter', 'latex');

    grid on;
    set(legend, 'Location', 'NorthEast');

    set(gca,'linewidth', axlabel_lindwidth)
    set(gca,'fontsize', axlabel_fontsize)
    xlim([1,num_mats]);
    xticks([1 20:20:100 num_mats]);
    
    ylim([1e-18 1])
    yticks(10.^(-18:3:0))

    % performance profile
    figure(2*j+2)

    ratio_max = 12;  
    markers = {'-v', '-o', '-s', '-^'};
    err_matix = [err_phifunc_j' err_phipade_opt_j' err_phipade_dft_j' err_expmblktri_j'];
    methods_names = {'\texttt{phi\_funm}', '\texttt{phipade\_opt}',...
        '\texttt{phipade\_dft}', '\texttt{expm\_blktri}'};

    myperfprof(err_matix, ratio_max, flip(mycolors(2:end,:),1), markers, lg_markersize, lg_lindwidth);
    legend(methods_names, 'interpreter', 'latex', 'FontSize', lg_fontsize);
    grid on;

    set(legend, 'Location', 'southeast');
    set(gca,'linewidth', axlabel_lindwidth)
    set(gca,'fontsize', axlabel_fontsize)

end

% figure for cost ratio
figure(11)  
    
cost_ratio_dft = cost_phifunc ./ cost_phipade_dft;
plot(1:num_mats, cost_ratio_dft(perm), 's', 'Color', [0.25, 0.25, 0.25], ...
    'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize);

grid on;

set(gca,'linewidth',axlabel_lindwidth)
set(gca,'fontsize',axlabel_fontsize)
xlim([1,num_mats]);
xticks([1 20:20:100 num_mats]);
    
ylim([0 1])

% title(sprintf('ratio of cost: \\texttt{phi\\_funm}/\\texttt{phipade\\_dft}'), 'FontSize', title_fontsize, 'Interpreter', 'latex');

figure(12)  
    
cost_ratio_opt = cost_phifunc ./ cost_phipade_opt;
plot(1:num_mats, cost_ratio_opt(perm), 's', 'Color', [0.25, 0.25, 0.25], ...
    'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize);

grid on;

set(gca,'linewidth',axlabel_lindwidth)
set(gca,'fontsize',axlabel_fontsize)
xlim([1,num_mats]);
xticks([1 20:20:100 num_mats]);
    
ylim([0 1])
% title(sprintf('ratio of cost: \\texttt{phi\\_funm}/\\texttt{phipade\\_opt}'), 'FontSize', title_fontsize, 'Interpreter', 'latex');

% figure for time ratio
figure(13)  
    
time_ratio_dft = time_phifunc ./ time_phipade_dft;
plot(1:num_mats, time_ratio_dft(perm), 's', 'Color', [0.25, 0.25, 0.25], ...
    'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize);

grid on;

set(gca,'linewidth',axlabel_lindwidth)
set(gca,'fontsize',axlabel_fontsize)
xlim([1,num_mats]);
xticks([1 20:20:100 num_mats]);
    
% ylim([0 1])

% title(sprintf('ratio of time: \\texttt{phi\\_funm}/\\texttt{phipade\\_dft}'), 'FontSize', title_fontsize, 'Interpreter', 'latex');

figure(14)  
    
time_ratio_opt = time_phifunc ./ time_phipade_opt;
plot(1:num_mats, time_ratio_opt(perm), 's', 'Color', [0.25, 0.25, 0.25], ...
    'LineWidth', lg_lindwidth, 'MarkerSize', lg_markersize);

grid on;

set(gca,'linewidth',axlabel_lindwidth)
set(gca,'fontsize',axlabel_fontsize)
xlim([1,num_mats]);
xticks([1 20:20:100 num_mats]);
    
% ylim([0 1])
% title(sprintf('ratio of time: \\texttt{phi\\_funm}/\\texttt{phipade\\_opt}'), 'FontSize', title_fontsize, 'Interpreter', 'latex');
