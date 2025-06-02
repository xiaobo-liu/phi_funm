function X = phi_func_ex(A, varargin)
% PHI_FUNC_EX Compute phi functions of a matrix for given indices using
% EXPM_MP in high precision arithmetic.
%
%   X = PHI_FUNC_EX(A, idx1, idx2, ..., idxN, d) computes the phi functions
%   of the matrix A for the indices specified in idx1, idx2, ..., idxN. The
%   last input, d, specifies the precision (number of digits) for the matrix
%   exponential computation. The output X is a cell array containing the
%   matrices corresponding to each index.
%
%   Inputs:
%       A       - A square matrix (n x n).
%       idx     - A vector or list of non-negative integers specifying the
%                 indices of the phi functions to compute. For example,
%                 idx = [1, 0, 3] computes phi_1(A), phi_0(A), and phi_3(A).
%       d       - (Optional) The number of digits of precision for the matrix
%                 exponential computation. Default is 200 if not provided.
%
%   Outputs:
%       X       - A cell array of size 1 x numel(indices), where X{i} contains
%                 the matrix corresponding to phi_j(A) for j = indices(i).
%
%   Example:
%       A = [1, 2; 3, 4];
%       idx = [1, 0, 3, 2]; % Compute phi_1(A), phi_0(A), phi_3(A), phi_2(A)
%       d = 100; % Precision of 100 digits
%       X = phi_func_ex(A, idx, d);
%
%       % Access results:
%       phi1 = X{1}; % phi_1(A)
%       phi0 = X{2}; % phi_0(A)
%       phi3 = X{3}; % phi_3(A)
%       phi2 = X{4}; % phi_2(A)
%
%   Notes:
%       1. The function uses a matrix exponential computation (expm_mp) to
%          evaluate the phi functions. Ensure that expm_mp is available in
%          your MATLAB environment.
%       2. If idx contains 0, the function computes phi_0(A) = expm(A).
%       3. The function supports both vector and scalar inputs for indices.
%
%   See also EXPM, KRON, SPEYE.

% Validate and process inputs
if isempty(varargin)
    error('At least one index must be provided.');
end

% Extract precision (d) if provided, otherwise use default value
if isscalar(varargin{end}) && mod(varargin{end},1)==0 && isfinite(varargin{end})...
        && varargin{end}>30 % d>30 needed, if provided
    d = varargin{end}; % Precision is the last input
    indices = [varargin{1:end-1}]; % All inputs except the last are indices
else
    d = 200; % Default precision
    indices = [varargin{:}]; % All inputs are indices
end

% Validate indices
if isempty(indices)
    error('At least one index must be provided.');
end

if ~(isnumeric(indices) && all(isfinite(indices)) && all(indices>= 0)...
        && all(mod(indices, 1) == 0))
    error('Indices must only contain nonnegative integers.')
end

p = max(indices);
n = size(A, 1);

% Compute phi functions
if p>0
    I = speye(n);
    E = zeros(n, n*p);
    E(:, 1:n) = I;
    J = tril(triu(ones(p), 1), 1);
    J = kron(J, speye(n));
    W = [A, E; zeros(n*p, n), J];
else % p=0
    W = A;
end

W = expm_mp(W, 'precision', d);

% Output format
X = cell(1, numel(indices)); % Store all matrices in a cell array
for i = 1:numel(indices)
    j = indices(i);
    X{i} = W(1:n, j*n+1:(j+1)*n);
end