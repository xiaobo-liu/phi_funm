% Test execution time against phipade.

% Create directories to store the results, if not exist
if ~exist('figs', 'dir'), mkdir('figs'); end
if ~exist('data', 'dir'), mkdir('data'); end

addpath('data','replication','figs')

rng default

pp = [1 4 7 10]; % here for p > 0  
pmax = max(pp);

num_pp = length(pp);

nn = [200 500 2500];
num_nn = length(nn);

u = eps('double')/2;
ids_min = 1;
[~, ids_max] = testmats_time();
num_mats = ids_max - ids_min + 1;

cost_phifunc = zeros(num_mats,1);
cost_phipade_opt = zeros(num_mats,1);
cost_phipade_dft = zeros(num_mats,1);

time_phifunc = zeros(num_mats,1);
time_phipade_opt = zeros(num_mats,1);
time_phipade_dft = zeros(num_mats,1);
        

for i=1:num_nn
    main_loop = tic; % record the time consumption for each n

    for k = ids_min:ids_max
        n = nn(i); % default matrix size
        fprintf('Running the test...Matrix id: %d\n', k);
        A = testmats_time(k,n);
        fprintf('Method      Time (seconds)\n');
        fprintf('---------------------------\n');

        I = eye(n);
        tic;
        [X_phifunc, ~, ~, cost_phifunc(k)] = phi_funm(A, pp); % X_phifunc{end} for p=0
        time_phifunc(k) = toc;
        fprintf('phi_funm     %10.6f\n', time_phifunc(k));

        % phipade degree and scaling optimized for cost 
        [deg_phipade, cost_phipade_opt(k)] = select_deg_phipade(A, pmax);
        tic;
        [G1, G2, G3, G4, G5, G6, G7, G8, G9, G10] = phipade(A, pmax, deg_phipade); 
        time_phipade_opt(k) = toc;
        fprintf('phipade_opt  %10.6f\n', time_phipade_opt(k));    

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
        fprintf('phipade_dft  %10.6f\n', time_phipade_dft(k)); 

        cost_phipade_dft(k) = phipade_default_cost(A, pmax);
    end
    cost_ratio_dft = cost_phifunc ./ cost_phipade_dft;
    cost_ratio_opt = cost_phifunc ./ cost_phipade_opt;
    time_ratio_dft = time_phifunc ./ time_phipade_dft;
    time_ratio_opt = time_phifunc ./ time_phipade_opt;
    [~, perm_dft] = sort(time_ratio_dft, 'descend');
    [~, perm_opt] = sort(time_ratio_opt, 'descend');

    % save the data for different n
    dataname = sprintf('data/time_%d.mat', n);    
    save(dataname, 'n', 'num_mats', 'cost_ratio_dft', 'cost_ratio_opt',...
        'time_ratio_dft', 'time_ratio_opt', 'perm_dft', 'perm_opt');

    fprintf('For n = %d, producing the results took %.2f minutes.\n',...
        n, toc(main_loop)/60);
end

    
%% load the data and plot

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

%% figure for time ratio phi_funm / phipade_dft

% sort the data according to the time ratio with the largest tested n
dataname = sprintf('data/time_%d.mat', max(nn));
load(dataname);
perm_dft1 = perm_dft;
perm_opt1 = perm_opt;

for i=num_nn:-1:1
    n = nn(i);
    dataname = sprintf('data/time_%d.mat', n);
    load(dataname);
    plot(1:num_mats, log2(1./time_ratio_dft(perm_dft1)), '-s', 'LineWidth', lg_lindwidth, ...
        'MarkerSize', lg_markersize);
    hold on
end
plot(1:num_mats, zeros(1,num_mats), '--', 'LineWidth', 2.5);
hold off

timelegend1 = sprintf('$n=%d$', nn(1));
timelegend2 = sprintf('$n=%d$', nn(2));
timelegend3 = sprintf('$n=%d$', nn(3));

legend(timelegend3, timelegend2, timelegend1, 'interpreter', 'latex', ...
    'Location', 'NW', 'FontSize', lg_fontsize);
set(gca,'linewidth',axlabel_lindwidth)
set(gca,'fontsize',axlabel_fontsize)
xlim([1,num_mats]);
xticks([1 10:10:50 num_mats]);
% title(sprintf('ratio of time: \\texttt{phi\\_funm}/\\texttt{phipade\\_dft}'), ...
        % 'FontSize', title_fontsize, 'Interpreter', 'latex');
grid on;

ylabel('$\log_2$ scale', 'interpreter', 'latex', 'FontSize', 1.5*lg_fontsize);
       
% ylim([0 1])
% yticks([0 0.2 0.4 0.8 1 1.25 2.5]);

figname = sprintf('figs/time_dft.eps');
exportgraphics(gca, figname, 'ContentType', 'vector');

%% figure for time ratio phi_funm / phipade_opt

figure(2)

for i=num_nn:-1:1
    n = nn(i);
    dataname = sprintf('data/time_%d.mat', n);
    load(dataname);
    plot(1:num_mats, log2(1./time_ratio_opt(perm_opt1)), '-s', 'LineWidth', lg_lindwidth, ...
        'MarkerSize', lg_markersize);
    hold on
end
plot(1:num_mats, zeros(1,num_mats), '--', 'LineWidth', 2.5);
hold off

timelegend1 = sprintf('$n=%d$', nn(1));
timelegend2 = sprintf('$n=%d$', nn(2));
timelegend3 = sprintf('$n=%d$', nn(3));

legend(timelegend3, timelegend2, timelegend1, 'interpreter', 'latex', ...
    'Location', 'NW', 'FontSize', lg_fontsize);
set(gca,'linewidth',axlabel_lindwidth)
set(gca,'fontsize',axlabel_fontsize)
xlim([1,num_mats]);
xticks([1 10:10:50 num_mats]);
% title(sprintf('ratio of time: \\texttt{phi\\_funm}/\\texttt{phipade\\_opt}'), ...
        % 'FontSize', title_fontsize, 'Interpreter', 'latex');
grid on;
       
% ylim([0 1])
% yticks([0 0.25 0.5 0.8 1 1.25 2]);

figname = sprintf('figs/time_opt.eps');
exportgraphics(gca, figname, 'ContentType', 'vector');