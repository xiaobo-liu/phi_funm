function [A, n_mats] = testmats(k, n)
%TESTMATS 108 nonnormal matrices from ANYMATRIX matrix collection, 
% including those from the matrices-expm collection at 
% https://github.com/xiaobo-liu/matrices-expm.

    n_mats = 108;
    if nargin < 1
        A = n_mats; % total number of matrices
        return;
    end
    if nargin == 2
        n = 20;
    end

    switch k
    case 01, sA='core/creation'; A=anymatrix(sA,n); % nilpotent
    case 02, sA='core/edelman27'; A=anymatrix(sA); % 27 by 27
    case 03, sA='core/gfpp'; A=anymatrix(sA,n); 
    case 04, sA='core/hessfull01'; A=anymatrix(sA,n); % rank=n-1
    case 05, sA='core/hess_orth'; A=anymatrix(sA,n);
    case 06, sA='core/nilpot_triang'; A=anymatrix(sA,n); % sparse
    case 07, sA='core/rschur'; A=anymatrix(sA,n);
    case 08, sA='core/stoch_cesaro'; A=anymatrix(sA,n); 
    case 09, sA='core/stoch_compan'; A=anymatrix(sA,n);
    case 10, sA='core/stoch_revtri'; A=anymatrix(sA,n);
    case 11, sA='core/tournament'; A=anymatrix(sA,n);
    case 12, sA='core/triminsval01'; A=anymatrix(sA,n);
    case 13, sA='core/vand'; A=anymatrix(sA,n);
    case 14, sA='gallery/binomial'; A=anymatrix(sA,n); 
    case 15, sA='gallery/chebspec'; A=anymatrix(sA,n);
    case 16, sA='gallery/chebvand'; A=anymatrix(sA,n);
    case 17, sA='gallery/chow'; A=anymatrix(sA,n);
    case 18, sA='gallery/circul'; A=anymatrix(sA,n);
    case 19, sA='gallery/clement'; A=anymatrix(sA,n);
    case 20, sA='gallery/condex'; A=anymatrix(sA,4,1); % 4 by 4
    case 21, sA='gallery/condex'; A=anymatrix(sA,3,2); % 3 by 3
    case 22, sA='gallery/condex'; A=anymatrix(sA,n,3); % lower triangular
    case 23, sA='gallery/cycol'; A=anymatrix(sA,n);
    case 24, sA='gallery/dorr'; A=anymatrix(sA,n); % sparse
    case 25, sA='gallery/dramadah'; A=anymatrix(sA,n);
    case 26, sA='gallery/dramadah'; A=anymatrix(sA,n,2);
    case 27, sA='gallery/dramadah'; A=anymatrix(sA,n,3);
    case 28, sA='gallery/forsythe'; A=anymatrix(sA,n);
    case 29, sA='gallery/forsythe'; A=anymatrix(sA,10,1e-10,1); 
    case 30, sA='gallery/frank'; A=anymatrix(sA,n);
    case 31, sA='gallery/gearmat'; A=anymatrix(sA,n);
    case 32, sA='gallery/grcar'; A=anymatrix(sA,n);
    case 33, sA='gallery/invhess'; A=anymatrix(sA,n);
    case 34, sA='gallery/invol'; A = anymatrix(sA,8)*8*pi; % \cite[exmp.~3]{hism03}
    case 35, sA='gallery/invol'; A = anymatrix(sA,13); A = triu(schur(A,'complex'),1); % nilpotent, triangular
    case 36, sA='gallery/jordbloc'; A=anymatrix(sA,n);
    case 37, sA='gallery/kahan'; A=anymatrix(sA,n);
    case 38, sA='gallery/leslie'; A=anymatrix(sA,n);
    case 39, sA='gallery/lesp'; A=anymatrix(sA,n);
    case 40, sA='gallery/lotkin'; A=anymatrix(sA,n);
    case 41, sA='gallery/neumann'; A=anymatrix(sA,ceil(sqrt(n))^2); % sparse, 25 by 25 for n=20;
    case 42, sA='gallery/normaldata'; A=anymatrix(sA,n,10);
    case 43, sA='gallery/orthog'; A=anymatrix(sA,n,-2);
    case 44, sA='gallery/parter'; A=anymatrix(sA,n);
    case 45, sA='gallery/randcolu'; A=anymatrix(sA,n);
    case 46, sA='gallery/randjorth'; A=anymatrix(sA,n);
    case 47, sA='gallery/rando'; A=anymatrix(sA,n,1);
    case 48, sA='gallery/rando'; A=anymatrix(sA,n,2);
    case 49, sA='gallery/rando'; A=anymatrix(sA,n,3);
    case 50, sA='gallery/rando'; A=10*triu(anymatrix(sA,n),1); % nilpotent, triangular.
    case 51, sA='gallery/randsvd'; A=anymatrix(sA,n,2);
    case 52, sA='gallery/randsvd'; A=anymatrix(sA,n,3);
    case 53, sA='gallery/randsvd'; A=anymatrix(sA,n,4);
    case 54, sA='gallery/randsvd'; A=anymatrix(sA,n,5);
    case 55, sA='gallery/randsvd'; A=anymatrix(sA,8,1e14); % \cite{fahi19}
    case 56, sA='gallery/redheff'; A=anymatrix(sA,n);
    case 57, sA='gallery/riemann'; A=anymatrix(sA,n);
    case 58, sA='gallery/sampling'; A=anymatrix(sA,n);
    case 59, sA='gallery/smoke'; A=anymatrix(sA,n);
    case 60, sA='gallery/smoke'; A=anymatrix(sA,n,1);
    case 61, sA='gallery/toeppen'; A=anymatrix(sA,n); % sparse
    case 62, sA='gallery/triw'; A=anymatrix(sA,n,-1);
    case 63, sA='gallery/triw'; A=anymatrix(sA,n,-2);
    case 64, sA='gallery/uniformdata'; A=anymatrix(sA,n,1e3);
    case 65, sA='gallery/wilk'; A=anymatrix(sA,3); % 3 by 3
    case 66, sA='gallery/wilk'; A=anymatrix(sA,4); % 4 by 3
    case 67, sA='matlab/pascal'; A=anymatrix(sA,n,1);
    case 68, sA='matlab/pascal'; A=anymatrix(sA,n,2);
    case 69, sA='nessie/spl0708b'; A=anymatrix(sA); % sparse, 41 by 41
       
    % matrices from matrices-expm collection
    case 70, sA='expm/alhi09r1'; A=anymatrix(sA);
    case 71, sA='expm/alhi09r2'; A=anymatrix(sA);
    case 72, sA='expm/alhi09r3'; A=anymatrix(sA);
    case 73, sA='expm/alhi09r4'; A=anymatrix(sA);
    case 74, sA='expm/dahi03'; A=anymatrix(sA);
    case 75, sA='expm/dipa00'; A=anymatrix(sA);
    case 76, sA='expm/edst04'; A=anymatrix(sA);
    case 77, sA='expm/eigt7'; A=anymatrix(sA);
    case 78, sA='expm/fahi19r1'; A=anymatrix(sA);
    case 79, sA='expm/fahi19r2'; A=anymatrix(sA);
    case 80, sA='expm/fahi19r3'; A=anymatrix(sA);
    case 81, sA='expm/fahi19r4'; A=anymatrix(sA);
    case 82, sA='expm/fasi7'; A=anymatrix(sA);
    case 83, sA='expm/jemc05r1'; A=anymatrix(sA);
    case 84, sA='expm/jemc05r2'; A=anymatrix(sA);
    case 85, sA='expm/kase99'; A=anymatrix(sA);
    case 86, sA='expm/kela89r1'; A=anymatrix(sA);
    case 87, sA='expm/kela89r2'; A=anymatrix(sA);
    case 88, sA='expm/kela98r1'; A=anymatrix(sA);
    case 89, sA='expm/kela98r2'; A=anymatrix(sA);
    case 90, sA='expm/kela98r3'; A=anymatrix(sA);
    case 91, sA='expm/kuda10'; A=anymatrix(sA);
    case 92, sA='expm/lara17r2'; A=anymatrix(sA);
    case 93, sA='expm/lara17r3'; A=anymatrix(sA);
    case 94, sA='expm/lara17r4'; A=anymatrix(sA);
    case 95, sA='expm/lara17r5'; A=anymatrix(sA);
    case 96, sA='expm/lara17r6'; A=anymatrix(sA);
    case 97, sA='expm/mopa03r1'; A=anymatrix(sA);
    case 98, sA='expm/mopa03r2'; A=anymatrix(sA);
    case 99, sA='expm/naha95'; A=anymatrix(sA);
    case 100, sA='expm/nies19'; A=anymatrix(sA);
    case 101, sA='expm/pang85r1'; A=anymatrix(sA);
    case 102, sA='expm/pang85r2'; A=anymatrix(sA);
    case 103, sA='expm/pang85r3'; A=anymatrix(sA);
    case 104, sA='expm/trem05'; A=anymatrix(sA);
    case 105, sA='expm/tsin13'; A=anymatrix(sA);
    case 106, sA='expm/ward77r1'; A=anymatrix(sA);
    case 107, sA='expm/ward77r3'; A=anymatrix(sA);
    case 108, sA='expm/ward77r4'; A=anymatrix(sA);   
    end
    A = full(A);
end