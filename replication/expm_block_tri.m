function [Fa,Fb,L] = expm_block_tri(A,B,E)

% EXPM_BLOCK_TRI Exponential of block triangular matrix.
% [Fa,Fb,L] = expm_block_tri(A,B,E) computes Fa = e^A, Fb = e^B, and
% L = D_{exp}(A,B,E), the (1,2) block of exp([A E;O B]), without explicitly
% computing the exponential of the block triangular matrix, exp([A E;O B]).
%
%   References:
%   A. H. Al-Mohy, Generalizing the Fr\'echet Derivative Algorithm for 
%      the Matrix Exponential, Submitted.
%
%   A. H. Al-Mohy and N. J. Higham, Computing the Fr\'echet Derivative of 
%      the Matrix Exponential, with an Application to Condition Number 
%      Estimation, SIAM J. Matrix Anal. Appl., 30(4), (2009), pp. 1639-1657.
%
%   A. H. Al-Mohy and N. J. Higham, A new scaling and squaring algorithm
%      for the matrix exponential, SIAM J. Matrix Anal. Appl., 31(3),
%      (2009), pp. 970-989.
%
%   Note: qtri_struct used by expm, a MATLAB internal function,
%         ...\toolbox\matlab\matfun\private\qtri_struct
%
%   Awad H. Al-Mohy, Sep. 3, 2024
     m_vals = [3 5 7 9 13];
             % ell_m for m=1:13.
     ell =   [%2.107342425540526e-008
              %3.555847927942764e-004
               1.081338577784837e-002  % m_vals = 3
              %6.486276265575182e-002
               1.998063206978949e-001  % m_vals = 5
              %4.373149570114990e-001
               7.834608472962045e-001  % m_vals = 7
              %1.234629393622150e+000
               1.782448623969279e+000  % m_vals = 9
              %2.416804616122239e+000
              %3.127459108657023e+000
              %3.904829673913844e+000
               4.740307543766806e+000];  % m_vals = 13

sizeA = length(A); sizeB = length(B);
classA = class(A); classB = class(B);
Ia = eye(sizeA,classA); Ib = eye(sizeB,classB); 
isrealA = isreal(A); isrealB  = isreal(B); isrealE = isreal(E);

eta = max(norm(A,1),norm(B,1));

recomputeDiagsAlow = false;
recomputeDiagsBlow = false;
schur_a = false; schur_b = false;

recomputeDiagsAup  = matlab.internal.math.isschur(A);
if ~recomputeDiagsAup
    recomputeDiagsAlow = matlab.internal.math.isschur(A.');
end
recomputeDiagsBup  = matlab.internal.math.isschur(B);
if ~recomputeDiagsBup
    recomputeDiagsBlow = matlab.internal.math.isschur(B.');
end
bigscalparam = (eta/ell(end) > 2^10);
if ~xor(recomputeDiagsAup , recomputeDiagsAlow) && bigscalparam
    schur_a = true; 
end
if ~xor(recomputeDiagsBup , recomputeDiagsBlow) && bigscalparam
    schur_b = true; 
end

if  schur_a
    [Qa,A] = schur(A,'complex'); 
    E = Qa'*E;
    recomputeDiagsAup = true;
end
if schur_b
    [Qb,B] = schur(B,'complex'); 
    E = E*Qb; 
    recomputeDiagsBup = true;
end

eta = max(norm(A,1),norm(B,1));

if recomputeDiagsAup
    blockformatAup = qtri_struct(A);
end
if recomputeDiagsAlow
    blockformatAlow = qtri_struct(A.');
end
if recomputeDiagsBup
    blockformatBup = qtri_struct(B);
end
if recomputeDiagsBlow
    blockformatBlow = qtri_struct(B.');
end            
        
if eta <= ell(end)
    % no scaling and squaring is required.
    for i = 1:length(m_vals)
        if eta <= ell(i)
            [Fa,Fb,L] = PadeApproximantOfDegree(m_vals(i));
            break;
        end
    end
else
    [t, s] = log2(eta/ell(end));
    s = s - (t == 0.5);    % Adjust s if eta/ell(end) is a power of 2.
    A = A/2^s; B = B/2^s; E = E/2^s; % Scaling
    [Fa,Fb,L] = PadeApproximantOfDegree(m_vals(end)); 
    if recomputeDiagsAup
        Fa = recompute_block_diag(A, Fa, blockformatAup);
    end
    if recomputeDiagsBup
        Fb = recompute_block_diag(B, Fb, blockformatBup);
    end
    if recomputeDiagsAlow
        Fa = recompute_block_diag(A.', Fa.', blockformatAlow);
        Fa = Fa.';
    end
    if recomputeDiagsBlow
        Fb = recompute_block_diag(B.', Fb.', blockformatBlow);
        Fb = Fb.';
    end
    for k= 1:s
        L = Fa*L + L*Fb; 
        Fa = Fa*Fa; Fb = Fb*Fb;
        if recomputeDiagsAup
            A = 2*A;
            Fa = recompute_block_diag(A, Fa, blockformatAup);
        end
        if recomputeDiagsBup
            B = 2*B;
            Fb = recompute_block_diag(B, Fb, blockformatBup);
        end
        if recomputeDiagsAlow
            A = 2*A;
            Fa = recompute_block_diag(A.', Fa.', blockformatAlow);
            Fa = Fa.';
        end
        if recomputeDiagsBlow
            B = 2*B;
            Fb = recompute_block_diag(B.', Fb.', blockformatBlow);
            Fb = Fb.';
        end
    end    
end
if schur_a
    Fa = Qa*Fa*Qa'; L = Qa*L;
    if isrealA, Fa = real(Fa); end
end

if schur_b
    Fb = Qb*Fb*Qb'; L = L*Qb';
    if isrealB, Fb = real(Fb); end
end
if isrealA && isrealB && isrealE, L = real(L); end

%%%%Nested Functions%%%%
 
    function [Fa,Fb,L] = PadeApproximantOfDegree(m)
        % PADEAPPROXIMANTOFDEGREE  Pade approximant to exponentials and the
        %   operator D_{exp}.
        %     
        %   [Fa,Fb,L] = PADEAPPROXIMANTOFDEGREE(M) is the degree M diagonal
        %   Pade approximant to EXP(A), EXP(B), and D_{exp}(A,B,E), where 
        %   M = 3, 5, 7, 9 or 13. That is, Fa = r_m(A), Fb = r_m(B), and 
        %   L = D_{r_m}(A,B,E).
        
        c = getPadeCoefficients;

        % Evaluate Pade approximant.
        switch m

            case {3, 5, 7, 9}
                m2 = (m+1)/2;      
                Apowers = cell(m2,1); Bpowers = cell(m2,1);
                Apowers{1} = Ia; Bpowers{1} = Ib;                 
                Apowers{2} = A*A; Bpowers{2} = B*B; 
                
                for j = 3:m2
                    Apowers{j} = Apowers{j-1}*Apowers{2};  
                    Bpowers{j} = Bpowers{j-1}*Bpowers{2};
                end
                
                Ua = zeros(sizeA,classA); Va = zeros(sizeA,classA);            
                Ub = zeros(sizeB,classB); Vb = zeros(sizeB,classB);
                for j = m+1:-2:2
                    Ua = Ua + c(j)*Apowers{j/2};  
                    Ub = Ub + c(j)*Bpowers{j/2};                   
                end                
                
                for j = m:-2:1
                    Va = Va + c(j)*Apowers{(j+1)/2}; 
                    Vb = Vb + c(j)*Bpowers{(j+1)/2};
                end
                
                   % Pade for D_{r_m} operator.
                   m3 = (m-1)/2;
                   M = cell(m3,1);  
                   M{1} = A*E + E*B;
                   for j = 2:m3
                       M{j} = Apowers{j}*M{1} + M{j-1}*Bpowers{2};
                   end
                   Lu = zeros(size(E),class(E)); Lv = zeros(size(E),class(E));
                   for j = m+1:-2:4
                       Lu = Lu + c(j)*M{j/2-1};
                   end
                   Lu = A*Lu + E*Ub;
                   for j = m:-2:3
                       Lv = Lv + c(j)*M{(j-1)/2};
                   end
                   Ua = A*Ua; Ub = B*Ub;

            case 13

                A2 = A*A; A4 = A2*A2; A6 = A2*A4;
                B2 = B*B; B4 = B2*B2; B6 = B2*B4;

                M2 = A*E + E*B;
                M4 = A2*M2+M2*B2;
                M6 = A2*M4+M2*B4;             
                               
                W1a = c(14)*A6 + c(12)*A4 + c(10)*A2;
                W2a = c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*Ia;
                Wa = A6*W1a + W2a;
                Ua = A*Wa;
                Z1a = c(13)*A6 + c(11)*A4 + c(9)*A2;
                Z2a = c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*Ia;
                Va= A6*Z1a + Z2a;


                W1b = c(14)*B6 + c(12)*B4 + c(10)*B2;
                W2b = c(8)*B6 + c(6)*B4 + c(4)*B2 + c(2)*Ib;
                Wb = B6*W1b + W2b;
                Ub = B*Wb;
                Z1b = c(13)*B6 + c(11)*B4 + c(9)*B2;
                Z2b = c(7)*B6 + c(5)*B4 + c(3)*B2 + c(1)*Ib;
                Vb= B6*Z1b + Z2b;

                Lw1 = c(14)*M6 + c(12)*M4 + c(10)*M2;
                Lw2 = c(8)*M6 + c(6)*M4 + c(4)*M2;
                Lw = A6*Lw1 + M6*W1b + Lw2;
                Lu = A*Lw + E*Wb;
                Lz1 = c(13)*M6 + c(11)*M4 + c(9)*M2;
                Lz2 = c(7)*M6 + c(5)*M4 + c(3)*M2;
                Lv = A6*Lz1 + M6*Z1b + Lz2;
        end

           %Fb = (-Ub+Vb)\(Ub+Vb);
           Fb = matlab.internal.math.nowarn.mldivide(Vb-Ub,2.*Ub) + Ib;

         %  [XL,XU] = lu(-Ua+Va); % LU does not preserve structure and deliver
         %    wrong result! Use A = anymatrix('matlab/pascal',30,1);
          % Fa = XU\(XL\(Ua+Va)); 
            Fa = (-Ua+Va)\(Ua+Va);
          % L = XU\(XL\(Lu+Lv + (Lu-Lv)*Fb)); 
            L = (-Ua+Va)\(Lu+Lv + (Lu-Lv)*Fb);

        
        function c = getPadeCoefficients
            % GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
            %    C = GETPADECOEFFICIENTS returns coefficients of numerator
            %    of [M/M] Pade approximant, where M = 3,5,7,9,13.
            switch m
                case 3
                    c = [120, 60, 12, 1]; 
                case 5
                    c = [30240, 15120, 3360, 420, 30, 1];
                case 7
                    c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
                case 9
                    c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
                         2162160, 110880, 3960, 90, 1];
                case 13
                    c = [64764752532480000, 32382376266240000, 7771770303897600, ...
                         1187353796428800,  129060195264000,   10559470521600, ...
                         670442572800,      33522128640,       1323241920,...
                         40840800,          960960,            16380,  182,  1];
            end
        end
    end


    function F = recompute_block_diag(T, F, block_struct)
    % recompute_block_diag Recomputes block diagonal of F = expm(T).
     n = length(T);
     for j = 1:n-1
        switch block_struct(j)
            case 0                
                % Not the start of a block, move on.
                continue;
            case 1
            % Start of a 2x2 triangular block.
                t11 = T(j,j);
                t22 = T(j+1,j+1);

                ave = (t11+t22)/2; df  = abs(t11-t22)/2;
                if max(ave,df) < log(realmax)                    
                   % Formula fine unless it overflows.
                    x12 = T(j,j+1)*exp(ave)*sinch((t22-t11)/2);
                else
                   % Revert to formula that can suffer cancellation.
                    x12 = T(j,j+1)*(exp(t22)-exp(t11))/(t22-t11);
                end
                F(j,j) = exp(t11);
                F(j,j+1) = x12;
                F(j+1,j+1) = exp(t22);
            case 2            
                % Start of a 2x2 quasi-triangular (full) block.
                a = T(j,j); b = T(j,j+1);
                c = T(j+1,j); d = T(j+1,j+1);
                delta = sqrt((a-d)^2 + 4*b*c)/2;
                expad2 = exp((a+d)/2);
                coshdelta = cosh(delta);
                sinchdelta = sinch(delta);
                F(j,j)     = expad2 .* (coshdelta + (a-d)./2.*sinchdelta);
                F(j+1,j)   = expad2 .* c .* sinchdelta;
               F(j,j+1)   = expad2 .* b .* sinchdelta; 
               F(j+1,j+1) = expad2 .* (coshdelta + (d-a)./2.*sinchdelta);
        end
     end
     % If last diagonal entry is not in a block it will have been missed.
     if block_struct(end) == 0
         F(n, n) = exp(T(n, n));
     end
     end

     function y = sinch(x)
     %sinch Returns sinh(x)/x.
     if x == 0
         y = 1;
     else
         y = sinh(x)/x;
     end
     end
end

%--------------------------------
% Subfunctions
%--------------------------------

function structure = qtri_struct(T)
%QTRI_STRUCT Block structure of a quasitriangular matrix T.
%
%   Let T be an n x n upper quasitriangular matrix then
%   structure is a list of numbers, of length n-1, 
%   where structure(j) encodes the block type of the j:j+1,j:j+1
%   diagonal block as one of the following.
%
%   0 - Not the start of a block.
%   1 - Start of a 2x2 triangular block.
%   2 - Start of a 2x2 quasi-triangular (full) block.

%   Nicholas J. Higham and Samuel D. Relton
%   Copyright 2014-2020 The MathWorks, Inc.

n = size(T,1);
if n == 1
    structure = 0;
    return;
elseif n == 2
    if T(2,1) == 0
        structure = 1;
        return;
    else
        structure = 2;
        return;
    end
end

j = 1;

structure = zeros(n-1, 1);

while j < n-1
    if T(j+1,j) ~= 0
        % Start of a 2x2 full block.
        structure(j:j+1) = [2; 0];
        j = j + 2; % Skip to next possible block start.
        continue;
    elseif T(j+1, j) == 0 && T(j+2, j+1) == 0
        % Start of a 2x2 triangular block.
        structure(j) = 1;
        j = j + 1;
        continue;
    else
        % Next block must start a 2x2 full block.
        structure(j) = 0;
        j = j + 1;
    end
end

% The n-1 entry has not yet been checked.
if T(n,n-1) ~= 0
    % 2x2 full block at the end.
    structure(n-1) = 2;
elseif (structure(n-2) == 0 || structure(n-2) == 1)
    structure(n-1) = 1;
end
end
%--------------------------------