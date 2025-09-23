%% This script runs all the tests.

anymatrix scan

format compact 
warning off 

% Create directories to store the results, if not exist
if ~exist('figs', 'dir'), mkdir('figs'); end
if ~exist('data', 'dir'), mkdir('data'); end

addpath('data','replication','figs')

global_time = tic();

%% SECTION - Test accuracy and cost against phipade and other algorithms.

fprintf('Running test_main...\n');
local_time = tic();
test_main
fprintf('done\t[%5.2f min]\n', toc(local_time) / 60);

close all 

%% SECTION - Test execution time against phipade.

fprintf('Running test_time...\n');
local_time = tic();
test_time
fprintf('done\t[%5.2f min]\n', toc(local_time) / 60);

close all

%% SECTION - Code profiling of phi_funm.

fprintf('Running test_profile...\n');
local_time = tic();
test_profile
fprintf('done\t[%5.2f min]\n', toc(local_time) / 60);

%% SECTION - Test on Hessenberg matrices.

fprintf('Running test_hess...\n');
local_time = tic();
test_hess
fprintf('done\t[%5.2f min]\n', toc(local_time) / 60);

%% End of test.

fprintf('All tests done\t[%5.2f min]\n', toc(global_time) / 60);
% exit