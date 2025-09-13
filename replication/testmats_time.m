function [A, n_mats] = testmats_time(k, n)
%TESTMATS_TIME 55 nonnormal matrices of variable size from ANYMATRIX matrix 
% collection.

    n_mats = 55;
    if nargin < 1
        A = n_mats; % total number of matrices
        return;
    end

switch k
    case 01, sA='core/creation'; A=anymatrix(sA,n); % nilpotent
    case 02, sA='core/gfpp'; A=anymatrix(sA,n); 
    case 03, sA='core/hessfull01'; A=anymatrix(sA,n); % rank=n-1
    case 04, sA='core/hess_orth'; A=anymatrix(sA,n);
    case 05, sA='core/nilpot_triang'; A=anymatrix(sA,n); % sparse
    case 06, sA='core/rschur'; A=anymatrix(sA,n);
    case 07, sA='core/stoch_cesaro'; A=anymatrix(sA,n); 
    case 08, sA='core/stoch_compan'; A=anymatrix(sA,n);
    case 09, sA='core/stoch_revtri'; A=anymatrix(sA,n);
    case 10, sA='core/tournament'; A=anymatrix(sA,n);
    case 11, sA='core/triminsval01'; A=anymatrix(sA,n);
    case 12, sA='core/vand'; A=anymatrix(sA,n);
    case 13, sA='gallery/chebspec'; A=anymatrix(sA,n);
    case 14, sA='gallery/chebvand'; A=anymatrix(sA,n);
    case 15, sA='gallery/chow'; A=anymatrix(sA,n);
    case 16, sA='gallery/circul'; A=anymatrix(sA,n);
    case 17, sA='gallery/clement'; A=anymatrix(sA,n);
    case 18, sA='gallery/condex'; A=anymatrix(sA,n,3); % lower triangular
    case 19, sA='gallery/cycol'; A=anymatrix(sA,n);
    case 20, sA='gallery/dorr'; A=anymatrix(sA,n); % sparse
    case 21, sA='gallery/dramadah'; A=anymatrix(sA,n);
    case 22, sA='gallery/dramadah'; A=anymatrix(sA,n,2);
    case 23, sA='gallery/dramadah'; A=anymatrix(sA,n,3);
    case 24, sA='gallery/forsythe'; A=anymatrix(sA,n);
    case 25, sA='gallery/frank'; A=anymatrix(sA,n);
    case 26, sA='gallery/gearmat'; A=anymatrix(sA,n);
    case 27, sA='gallery/grcar'; A=anymatrix(sA,n);
    case 28, sA='gallery/invhess'; A=anymatrix(sA,n);
    case 29, sA='gallery/jordbloc'; A=anymatrix(sA,n);
    case 30, sA='gallery/kahan'; A=anymatrix(sA,n);
    case 31, sA='gallery/leslie'; A=anymatrix(sA,n);
    case 32, sA='gallery/lesp'; A=anymatrix(sA,n);
    case 33, sA='gallery/lotkin'; A=anymatrix(sA,n);
    case 34, sA='gallery/normaldata'; A=anymatrix(sA,n,10);
    case 35, sA='gallery/orthog'; A=anymatrix(sA,n,-2);
    case 36, sA='gallery/parter'; A=anymatrix(sA,n);
    case 37, sA='gallery/randcolu'; A=anymatrix(sA,n);
    case 38, sA='gallery/randjorth'; A=anymatrix(sA,n);
    case 39, sA='gallery/rando'; A=anymatrix(sA,n,1);
    case 40, sA='gallery/rando'; A=anymatrix(sA,n,2);
    case 41, sA='gallery/rando'; A=anymatrix(sA,n,3);
    case 42, sA='gallery/rando'; A=10*triu(anymatrix(sA,n),1); % nilpotent, triangular.
    case 43, sA='gallery/randsvd'; A=anymatrix(sA,n,2);
    case 44, sA='gallery/randsvd'; A=anymatrix(sA,n,3);
    case 45, sA='gallery/randsvd'; A=anymatrix(sA,n,4);
    case 46, sA='gallery/randsvd'; A=anymatrix(sA,n,5);
    case 47, sA='gallery/redheff'; A=anymatrix(sA,n);
    case 48, sA='gallery/riemann'; A=anymatrix(sA,n);
    case 49, sA='gallery/sampling'; A=anymatrix(sA,n);
    case 50, sA='gallery/smoke'; A=anymatrix(sA,n);
    case 51, sA='gallery/smoke'; A=anymatrix(sA,n,1);
    case 52, sA='gallery/toeppen'; A=anymatrix(sA,n); % sparse
    case 53, sA='gallery/triw'; A=anymatrix(sA,n,-1);
    case 54, sA='gallery/triw'; A=anymatrix(sA,n,-2);
    case 55, sA='gallery/uniformdata'; A=anymatrix(sA,n,1e3);
    % case 13, sA='gallery/binomial'; A=anymatrix(sA,n); % overflow issue
    % case 56, sA='matlab/pascal'; A=anymatrix(sA,n,1); % NaN for phipade
    % case 58, sA='matlab/pascal'; A=anymatrix(sA,n,2); % overflow issue
      
    A = full(A);
end