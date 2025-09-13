function [X,s,m,cost] = phi_funm(A,varargin)
%PHI_FUNM Compute phi functions of a matrix for given indices.
%
%   X = PHI_FUNM(A, idx1, idx2, ..., idxN) computes the phi 
%   functions of the matrix A for the indices specified in 
%                      idx1, idx2, ..., idxN. 
%   The output X is a cell array containing the matrices corresponding to 
%   each index.
%
%   Inputs:
%       A      - A square matrix (n x n).
%       idx    - A vector or list of non-negative integers specifying the
%                indices of the phi functions to compute. For example,
%                idx = [1, 0, 3] computes phi_1(A), phi_0(A), and phi_3(A).
%
%   Outputs:
%       X      - A cell array of size 1 x numel(indices), where X{i} contains
%                the matrix corresponding to phi_j(A) for j = indices(i).
%       s      - The scaling parameter.
%       m      - The degree of the Pade approximant to phi_p.
%       cost   - The overall computational cost measured in the equivalent 
%                number of matrix products.
%
% Reference:
%
%    Awad H. Al-Mohy and Xiaobo Liu. Computing Matrix $\varphi$-Functions 
%    Arising in Exponential Integrators. ArXiv:2506.01193 [math.NA], June 2025.
%
% Version Date: June 1, 2025

indices = [varargin{:}];
p = max(indices); 
if p == 0
    error('Use the MATLAB built-in function expm for computing \varphi_0(A) = e^A.');
end
% Check if A is in upper or lower Schur form
recomputeDiagsExpUp  = matlab.internal.math.isschur(A);

if ~recomputeDiagsExpUp
    recomputeDiagsExpLow = matlab.internal.math.isschur(A.');
else
    recomputeDiagsExpLow = false;
end
% Capture the structure of  the 2-by-2 diagonal blocks 
if recomputeDiagsExpUp
    blockformatAup = qtri_struct(A);
end
if recomputeDiagsExpLow
    blockformatAlow = qtri_struct(A.');
end 

classA = class(A);
n = length(A); I = eye(n,classA);

[m,s,tau,cost] = select_parameters_phi(A,p);  

A = A./2^s;                               
[N, D] = pade_coef(m,p);                   
[Nm,Dm] = Paterson_Stockmeyer(A,N{p},D{p},tau); % Paterson-Stockmeyer method

Rm = zeros(n,n,p+1,classA);
% Pade approximant to phi_p
Rm(:,:,p+1) = matlab.internal.math.nowarn.mldivide(Dm,Nm); % Dm\Nm

coef = cumprod(1:p-1);  % Computes [1!, 2!, 3!, ..., (p-1)!]
coef = [1, coef];       % Add 0! = 1 at the beginning
for k=p:-1:1 
    Rm(:,:,k) = A * Rm(:,:,k+1) + I./coef(k);
end
% Recompute the 2-by-2 diagonal blocks of EXPM
if recomputeDiagsExpUp
   Rm(:,:,1) = recompute_block_diag(A, Rm(:,:,1), blockformatAup);
end    
if recomputeDiagsExpLow
    Rm(:,:,1) = recompute_block_diag(A.', Rm(:,:,1).', blockformatAlow);
    Rm(:,:,1) = Rm(:,:,1).';
end

if s>0, flipcoef = 1./flip(coef); end
for i=1:s 
    for j=p:-1:1 
        Rm(:,:,j+1) = (Rm(:,:,1)*Rm(:,:,j+1) + ...
            sum(bsxfun(@times, Rm(:,:,2:j+1), reshape(flipcoef(p-j+1:end), [1,1,j])),3))./2^j;
    end         
        Rm(:,:,1) = Rm(:,:,1)*Rm(:,:,1);   % EXPM     
    
    % Recompute the 2-by-2 diagonal blocks of EXPM
    if recomputeDiagsExpUp
        Rm(:,:,1) = recompute_block_diag(2^i.*A, Rm(:,:,1), blockformatAup);
    end
    if recomputeDiagsExpLow
        Rm(:,:,1) = recompute_block_diag(2^i.*A.', Rm(:,:,1).', blockformatAlow);
        Rm(:,:,1) = Rm(:,:,1).';
    end 
end      
% Output format
X = cell(1, length(indices));
for i = 1:length(indices)
    j = indices(i);
    X{i} = Rm(:,:,j+1);
end
end

%--------------------------------
% Subfunctions

%--------------------------------
function [m, s, tau, cost] = select_parameters_phi(A,p)

% SELECT_PARAMETERS_PHI Parameters for the phi-functions algorithm. 
%
%  [M,S,TAU] = SELECT_PARAMETERS_PHI(A,P) computes the optimal scaling 
%  parameter S, the degree of the [M/M] Pade approximant to phi_p, and the 
%  number of the powers of the matrix A, TAU. 
%
%  The number of matrix multiplications requied for evaluating a [m/m] Pade
%  approximant, N_m(A)/D_m(A), where N_m and D_m are polynomials of 
%  degree m, using Paterson–Stockmeyer scheme, is 
%
%      pi_m(tau) = tau - 1 + 2*floor(m/tau) - 2*( floor(m/tau)==m/tau ), 
%
%  where tau is the degree of required powers of A: A^2, A^3, ..., A^tau.
%  Theoretically, tau_* that minimizes pi_m(tau) is either 
%  tau = floor(sqrt(2*m)) or tau = ceil(sqrt(2*m)).
%  The optimal degrees for the scheme is gigen explicitly by
%
%              m_i = floor( (i + 3)^2/8 ), i = 0, 1, 2, ...
%
%  At these optimal degrees, pi_{m_i}(tau_*) = i, see ref. [2].
%
%  The total number of matrix multiplications needed for evaluating the 
%  phi-functions, phi_0, phi_1, ..., phi_p, using optimal orders, is
%
%             C(m_i,r) =: pi_{m_i}(tau_*) + p + 4/3 + max(s0,t)*(p+1),
%                      =        i         + p + 4/3 + max(s0,t)*(p+1)
%
%  where s0 is a scaling parameter such that 2^(-s0)*alpha_r <= theta(m_i,p)
%  and t is another scaling parameter related to the first term of the
%  backward error series.
%  This code computes m and s = max(s0,t) that minimize C(m_i,r).
%
% References:
%
% [1] Awad H. Al-Mohy and Xiaobo Liu. Computing Matrix $\varphi$-Functions 
%        Arising in Exponential Integrators. ArXiv:2506.01193 [math.NA], June 2025.
% [2] Massimiliano Fasi. Optimality of the Paterson–Stockmeyer method for 
%        evaluating matrix polynomials and rational matrix functions. 
%        Linear Algebra Appl., 574:182–200, 2019.

%load('theta_m_p.mat','theta');
theta = parameter_theta_mp();
m_max = 12; % Larger m (up to 20) increases cond(D_m(2^{-s}A)) and may reduce accuracy.

if p > 7, theta(:,p) = theta(:,7); end % For the stability of Pade approx.       

i_max = ceil(sqrt(8*(m_max + 1)) - 3) - 1;
m_max = floor( (i_max + 3)^2/8 );
phat = p;
if theta(m_max,p) < 1, phat = 0; end % In case a use changes m_max ( < 7 )

% The largest r such that 2*m_max + phat +1 >= r(r-1)
r_max = floor( ( 1 + sqrt(1 + 4*(2*m_max + phat + 1)) )/2 );

eta = zeros(r_max,1); alpha = zeros(r_max-1,1);
for j = 1:r_max
    c = normAm(A,j+1); % 1-norm estimation of norm(A^(j+1),1)
    c = c^(1/(j+1));
    eta(j) = c;
end
for j = 1:r_max-1
    alpha(j) = max(eta(j),eta(j+1));
end

Cost_matrix = zeros(i_max+1,r_max-1); % cost matrix
for i= 0:i_max 
    m_i = floor( (i + 3)^2/8 );
    phat = p;
    if theta(m_i,p) < 1, phat = 0; end
    t = ell(A,m_i,p,phat); % scaling paramenter via the first term.
    for r = 2:r_max
        if 2*m_i + phat + 1 >= r*(r - 1)
            Cost_matrix(i+1, r-1) = i + p + ...  % remove 4/3 to have integer elements
                max(ceil(log2(alpha(r - 1) / theta(m_i, p))), t) * (p + 1);           
        end
    end
end

pve_elmnts = Cost_matrix( Cost_matrix > 0 ); % All positive elements.
min_pve_elmnts = min(pve_elmnts);

cost = min_pve_elmnts + 4/3;

% The row and column indices of the smallest positive element
%[row, r_star] = find(Cost_matrix == min_pve_elmnts, 1, 'last');
[row, ~] = find(Cost_matrix == min_pve_elmnts, 1, 'last');
i_star = row - 1; % The index is shifted above by 1
% The optimal degree of Pade approx.
m = floor( (i_star + 3)^2/8 );

% The optimal scaling parameter
s = (min_pve_elmnts  - i_star - p)/(p+1);
% Number of required powers of A to evaluate Pade of phi_p.
tau = floor(sqrt(2*m));
if tau - 1 + 2*floor(m/tau) - 2*( floor(m/tau) == m/tau ) ~= i_star
    tau = ceil(sqrt(2*m));
end
end

%--------------------------------
function t = ell(T, m, p,phat)
%ell Function needed to correct the scaling paramter using the first term
% of the backward error function h_{m,p}.
normT = norm(T,1);
t0 = 0;     
if normT > 1, t0 = log2(normT); end
T = T/2^t0;         % prescale T to avoid overflow of alpha
normT = normT/2^t0;
delta = (p-1)*(p-phat)/p + 1;
coeff = ( factorial(m+p)/factorial(2*m+p) )*...
        ( factorial(m)/factorial(2*m+p+1) );

scaledT = (coeff/(normT^delta)).^(1/(2*m+p+1)) .* abs(T);
alpha = normAm(scaledT,2*m+p+1);
t = log2(2*alpha/eps(class(alpha)))/(2*m+p+1-delta) + t0; % undo prescaling
t = max( ceil(t) , 0);
end

%--------------------------------
function [N, D] = pade_coef(m, p)
% Re-normalised [m/m] Padé coefficients for the \phi_k functions.
%
% This code is based on the explicit [m/m] Padé formulae derived by W.
% Wright, but uses recurrence relations to reduce the computational
% complexity.  In particular, the coefficients n_i of the numerator
% polynomial N_m^k are defined by
%
%    n_i = \sum_{j=0}^i a_{ij} = sum(tril(A), 2)_i
%
% for i=0:m, in which
%
%    a_{ij} = (2m + k - j)! (-1)^j / (j! (m-j)! (k + i - j)!)
%           = -(m+1 - j) (k+1+i - j) / (j (2m+k+1 - j)) a_{i,j-1}
%
% for j=1:i.  Similar recurrence relations may be derived for the other
% coefficients, and for the denominator polynomial D_m^k.
%
% We note that roundoff errors unfortunately affects the accuracy of the
% coefficients.  However, as the errors are generally in the order of
% 1-5 ULP, we do not implement more accurate evaluation routines at this
% time.
%
%  Reference:
%
% H. Berland, B. Skaflestad, and W. M. Wright}, EXPINT--a MATLAB
%  package for exponential integrators, ACM Trans. Math. Softw., 33 (2007),
%  p.~4–es, https://doi.org/10.1145/1206040.1206044.

n1 = prod(m + 1 : 2*m + 1);     % (2m + 1)! / m!
d1 = n1;

i = 1:m;
[J, I] = meshgrid(i);   % MESHGRID gives wrong order for this purpose

N = cell([p, 1]);
D = cell([p, 1]);
A = zeros(m + 1);

k = 1;
while k <= p
   A(:, 1) = n1 .* cumprod([1, 1 ./ (k + i)]) .';
   A(2:end, 2:end) = - (m + 1 - J) .* (k + 1 + I - J) ./ ...
                       ((2*m + k + 1 - J) .* J);

   N{k} = sum(tril(cumprod(A, 2)), 2);
   D{k} = d1 .* cumprod([1, -(m + 1 - i) ./ (i .* (2*m + k + 1 - i))]) .';

   k = k + 1;

   n1 = n1 * (2*m + k) / k;
   d1 = d1 * (2*m + k);
end
end

%--------------------------------
function [N, D] = Paterson_Stockmeyer(A,Ncoeffs, Dcoeffs, tau)
    % PATERSON_STOCKMEYER Efficiently evaluates two matrix polynomials of
    % the same degree.
    %
    %   [N, D] = Paterson_Stockmeyer(Ncoeffs, Dcoeffs, A, tau)
    %   evaluates the matrix polynomials:
    %       N(A) = c0_1*I + c1_1*A + c2_1*A^2 + ... + cm_1*A^m
    %       D(A) = c0_2*I + c1_2*A + c2_2*A^2 + ... + cm_2*A^m
    %   using the Paterson-Stockmeyer method with block size tau.
    %
    % Inputs:
    %   Ncoeffs - Coefficients of N.
    %   Dcoeffs - Coefficients of D.
    %   A       - The matrix at which to evaluate the polynomials.
    %   tau     - The block size (typically sqrt(2m)).
    %
    % Outputs:
    %   N - The value of the first matrix polynomial N(A).
    %   D - The value of the second matrix polynomial D(A).

m = length(Ncoeffs) - 1; % Degree of the polynomials

% The number of blocks
nu = floor(m / tau);

% Precompute powers of A up to A^tau and store them in a 3D array
n = size(A,1);
A_powers = zeros(n, n, tau+1); % 3D array to store A^0, A^1, ..., A^tau
A_powers(:, :, 1) = eye(n); % A^0
A_powers(:, :, 2) = A;      % A^1
for i = 3:tau+1
    A_powers(:, :, i) = A_powers(:, :, i-1) * A; % A^2, ..., A^tau
end
Atau = A_powers(:, :, tau+1); % A^tau

% Initialize results for both polynomials
N = zeros(n);
D = zeros(n);

% Evaluate the polynomials in blocks
for i = 0:nu
% Determine the start and end indices of the current block
    start_idx = i * tau + 1;
    end_idx = min((i + 1) * tau, m + 1);
   % Extract coefficients for the current block
    blk_Ncoeffs = Ncoeffs(start_idx:end_idx);
    blk_Dcoeffs = Dcoeffs(start_idx:end_idx);

   % Extract corresponding powers of A from the 3D array
    block_powers = A_powers(:, :, 1:length(blk_Ncoeffs));

   % Evaluate the block polynomials B_i^{[p]}(A) for both polynomials
    blkvalN = sum(bsxfun(@times,block_powers , reshape(blk_Ncoeffs, 1, 1, [])), 3);                    
    blkvalD = sum(bsxfun(@times,block_powers , reshape(blk_Dcoeffs, 1, 1, [])), 3);   

    % Multiply by (A^tau)^i and add to the results
    if i == 0
       N = N + blkvalN; % (A^s)^0 = I
       D = D + blkvalD; % (A^s)^0 = I
    else
        N = N + blkvalN * Atau;
        D = D + blkvalD * Atau;
        Atau = Atau * A_powers(:, :, tau+1); % Update Atau to (A^tau)^(i+1)
    end
end
end

%--------------------------------
function F = recompute_block_diag(T, F, block_struct)
%RECOMPUTE_BLOCK Recomputes block diagonal of F = expm(T).
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

%--------------------------------
function y = sinch(x)
%sinch Returns sinh(x)/x.
if x == 0
   y = 1;
else
   y = sinh(x)/x;
end
end

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
function [c,mv] = normAm(A,m)
%NORMAM   Estimate of 1-norm of power of matrix.
%   NORMAM(A,m) estimates norm(A^m,1).
%   If A has nonnegative elements the estimate is exact.
%   [C,MV] = NORMAM(A,m) returns the estimate C and the number MV of
%   matrix-vector products computed involving A or A^*.
%
%   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
%   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
%   970-989, 2009.
%   Awad H. Al-Mohy and Nicholas J. Higham, September 7, 2010.

t = 1; % Number of columns used by NORMEST1.

n = length(A);
if isequal(A,abs(A))
    e = ones(n,1);
    for j=1:m         % for positive matrices only
        e = A'*e;
    end
    c = norm(e,inf);
    mv = m;
else
    [c,v,w,it] = normest1(@afun_power,t);
    mv = it(2)*t*m;
end

  function Z = afun_power(flag,X)
       %AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.

       if isequal(flag,'dim')
          Z = n;
       elseif isequal(flag,'real')
          Z = isreal(A);
       else

          [p,q] = size(X);
          if p ~= n, error('Dimension mismatch'), end

          if isequal(flag,'notransp')
             for i = 1:m, X = A*X; end
          elseif isequal(flag,'transp')
             for i = 1:m, X = A'*X; end
          end

          Z = X;

       end

  end
end

%-----------------------------------
function theta = parameter_theta_mp()
%PARAMETER_THETA_MP Returns a 20-by-10 matrix of theta_{m,p} values.
%
%   theta = PARAMETER_THETA_MP() returns a 20-by-10 matrix theta
%   containing precomputed parameter values \theta_{m,p} used for
%   scaling in phi-function approximations.
theta = [ ...
    % Column 1
     1.9994634524084072e-05 ...
     3.8062018282832713e-03 ...
     3.9716360056613338e-02 ...
     1.5442675548312682e-01 ...
     3.7980898016147974e-01 ...
     7.2617719570332095e-01 ...
     1.1898666361923196e+00 ...
     1.7605812331512907e+00 ...
     2.4258189585472332e+00 ...
     3.1731134567937489e+00 ...
     3.9910252013298151e+00 ...
     4.8694854897845783e+00 ...
     5.7998279045478851e+00 ...
     6.7746819964584413e+00 ...
     7.7878170515724348e+00 ...
     8.8339782148338468e+00 ...
     9.9087338231338435e+00 ...
     1.1008341035221061e+01 ...
     1.2129631147167878e+01 ...
     1.3269913405576162e+01; ...
    % Column 2
     3.7631213142601604e-05 ...
     6.0902062861257263e-03 ...
     5.8069688868806917e-02 ...
     2.1278117034577634e-01 ...
     5.0146249765870066e-01 ...
     9.2819101596462739e-01 ...
     1.4469917284172122e+00 ...
     2.0609071947420161e+00 ...
     2.7623053687512256e+00 ...
     3.5393251225873685e+00 ...
     4.3813599805888455e+00 ...
     5.2791998702489158e+00 ...
     6.2249713216199041e+00 ...
     7.2119972401875811e+00 ...
     8.2346349750856938e+00 ...
     9.2881195676692201e+00 ...
     1.0368423129623023e+01 ...
     1.1472133539821815e+01 ...
     1.2596352058386916e+01 ...
     1.3738607964511967e+01; ...
    % Column 3
     7.3660064160451633e-05 ...
     9.8696827467796150e-03 ...
     8.5342200767598173e-02 ...
     2.9371996708854947e-01 ...
     6.6210338244568179e-01 ...
     1.1591052927815408e+00 ...
     1.7187282676344413e+00 ...
     2.3714800303152570e+00 ...
     3.1051747305115041e+00 ...
     3.9086764034913664e+00 ...
     4.7722044322066060e+00 ...
     5.6873445011759571e+00 ...
     6.6469348065071401e+00 ...
     7.6449101286829073e+00 ...
     8.6761421967665289e+00 ...
     9.7362927928661431e+00 ...
     1.0821685340709736e+01 ...
     1.1929195701618417e+01 ...
     1.3056160721729096e+01 ...
     1.4200302290386791e+01; ...
    % Column 4
     1.4973317297025854e-04 ...
     1.6211831146383013e-02 ...
     1.2612151517169362e-01 ...
     4.0617647304246707e-01 ...
     8.7398587777443537e-01 ...
     1.4012982671152012e+00 ...
     2.0024009817133357e+00 ...
     2.6900789174368360e+00 ...
     3.4527057639100476e+00 ...
     4.2799119378223347e+00 ...
     5.1627038370305609e+00 ...
     6.0933928000926798e+00 ...
     7.0654517332116464e+00 ...
     8.0733547678005788e+00 ...
     9.1124243895036638e+00 ...
     1.0178695474892383e+01 ...
     1.1268798592552891e+01 ...
     1.2379861767254168e+01 ...
     1.3509428665638305e+01 ...
     1.4655390835644507e+01; ...
    % Column 5
     3.1524433337711818e-04 ...
     2.6984240563843326e-02 ...
     1.8736524835661617e-01 ...
     5.6238430023209962e-01 ...
     1.1108284429014510e+00 ...
     1.6570386251184455e+00 ...
     2.2957295231389230e+00 ...
     3.0148775983336011e+00 ...
     3.8035234182863555e+00 ...
     4.6520580625414274e+00 ...
     5.5522206158634058e+00 ...
     6.4969779362477960e+00 ...
     7.4803674046975299e+00 ...
     8.4973389493985039e+00 ...
     9.5436112198248839e+00 ...
     1.0615546677548995e+01 ...
     1.1710045790427412e+01 ...
     1.2824458619109002e+01 ...
     1.3956511465304834e+01 ...
     1.5104246216977073e+01; ...
    % Column 6
     6.8552099837644348e-04 ...
     4.5475796838180302e-02 ...
     2.7955116495245524e-01 ...
     7.7878835057542650e-01 ...
     1.3375805277521873e+00 ...
     1.9240266024163379e+00 ...
     2.5968019781751672e+00 ...
     3.3443898815994957e+00 ...
     4.1565381842964753e+00 ...
     5.0243666244242222e+00 ...
     5.9402869155388993e+00 ...
     6.8978550856394785e+00 ...
     7.8916102694998589e+00 ...
     8.9169229988464433e+00 ...
     9.9698613731505183e+00 ...
     1.1047076803237459e+01 ...
     1.2145708155762449e+01 ...
     1.3263302077220441e+01 ...
     1.4397747062985934e+01 ...
     1.5547218984875734e+01; ...
    % Column 7
     1.5357294906993542e-03 ...
     7.7492732598553313e-02 ...
     4.1822735418681667e-01 ...
     1.0464245027100287e+00 ...
     1.5770955744229305e+00 ...
     2.2002921672563200e+00 ...
     2.9040336511997444e+00 ...
     3.6774154946461342e+00 ...
     4.5108926866861268e+00 ...
     5.3962683334405188e+00 ...
     6.3265666055677787e+00 ...
     7.2958719120586348e+00 ...
     8.2991698538199703e+00 ...
     9.3322038606787086e+00 ...
     1.0391351525913848e+01 ...
     1.1473520338548882e+01 ...
     1.2576060822320501e+01 ...
     1.3696694607678852e+01 ...
     1.4833455011957760e+01 ...
     1.5984637960610213e+01; ...
    % Column 8
     3.5357368946407606e-03 ...
     1.3324895779231027e-01 ...
     6.2587028003519907e-01 ...
     1.2572921799132364e+00 ...
     1.8273848586096524e+00 ...
     2.4841706972729956e+00 ...
     3.2161214838471062e+00 ...
     4.0129905619037700e+00 ...
     4.8659163980540407e+00 ...
     5.7673348160254694e+00 ...
     6.7108251879721594e+00 ...
     7.6909456666904150e+00 ...
     8.7030799466790363e+00 ...
     9.7433031993539103e+00 ...
     1.0808268293369911e+01 ...
     1.1895110742359616e+01 ...
     1.3001369927621694e+01 ...
     1.4124924034990169e+01 ...
     1.5263936358537684e+01 ...
     1.6416810943368130e+01; ...
    % Column 9
     8.3457890280622341e-03 ...
     2.3046256043625210e-01 ...
     9.3359704439865621e-01 ...
     1.4803508614519840e+00 ...
     2.0867311457667101e+00 ...
     2.7742685942850143e+00 ...
     3.5319998551192429e+00 ...
     4.3503444488739316e+00 ...
     5.2210881567935115e+00 ...
     6.1372481669903998e+00 ...
     7.0929061730821310e+00 ...
     8.0830455254242484e+00 ...
     9.1034058521773762e+00 ...
     1.0150358556572542e+01 ...
     1.1220802409044003e+01 ...
     1.2312076916038997e+01 ...
     1.3421890788016032e+01 ...
     1.4548262963142152e+01 ...
     1.5689473955602649e+01 ...
     1.6844025647301478e+01; ...
    % Column 10
     2.0138828089282260e-02 ...
     3.9889911045491461e-01 ...
     1.1616793320890249e+00 ...
     1.7138710031853250e+00 ...
     2.3536689190936406e+00 ...
     3.0694262945895781e+00 ...
     3.8508005672569903e+00 ...
     4.6888632793636376e+00 ...
     5.5760056404345404e+00 ...
     6.5057767447444323e+00 ...
     7.4727126397309815e+00 ...
     8.4721790248909414e+00 ...
     9.5002347684515378e+00 ...
     1.0553516823928600e+01 ...
     1.1629144570605774e+01 ...
     1.2724640835984562e+01 ...
     1.3837866852904227e+01 ...
     1.4966968689718191e+01 ...
     1.6110333058994453e+01 ...
     1.7266550769620974e+01 ]';
end